#include "mex.h"
#include "class_handle.hpp"
#include <iostream>
#include <AdaGrad.h>



AdaGrad::AdaGrad(Matrix<double, num_state_variables, num_state_variables> A, Matrix<double, num_state_variables, num_manipulated_variables>  B,
	Matrix<double, num_state_variables, num_state_variables> Bd, Matrix<double, num_state_variables, num_state_variables> Q,
	Matrix<double, num_state_variables, num_state_variables> Q_final, Matrix<double, num_manipulated_variables, num_manipulated_variables> R,
	Matrix<double, num_manipulated_variables, num_manipulated_variables> R_delta,
	Matrix<double, num_state_variables, 1> disturbance, int num_params, int pred_horizon, int control_horizon,
	Matrix<double, num_manipulated_variables, num_manipulated_variables> scale_MV, Matrix<double, num_state_variables, num_state_variables> scale_OV) : num_params_(num_params), pred_horizon(pred_horizon), control_horizon(control_horizon), A_(A), B_(B), Bd_(Bd), Q_(Q), Q_final_(Q_final), R_(R), R_delta_(R_delta), disturbance_(disturbance)
{
	//this->insecure_ = this->Bd_ * disturbance;
	x_states.setZero();
	deriv_wrt_u.setZero();
	u.setZero();
	u_past.setZero();
	u_current.setZero();
	u_ss_.setZero();
	lambdas_x.setZero();
	lambdas_u.setZero();
	lambdas_u_ref.setZero();
	lambdas_u_ref.setZero();
	x0_.setZero();
	change_x.setZero();
	x_ss_.setZero();
	gradients.setZero();
	gradient_past.setZero();
	epsilon = 0.4;
	delta = 1e-7;
	residuum = 0.0;
	max_iter = 50;
	alfa = 200;
	saturation_count = 0;
	min_residuum = 1e-5;
	count_jacobians = 0;
	scale_MV_inv = scale_MV.inverse();
	scale_OV_inv = scale_OV.inverse();
	A_pow_B_cache.setZero();
	A_pow_B_cache.block(0, 0, A_.rows(), B_.cols()) = MatrixXd::Identity(A_.rows(), A_.cols())* B_;

	for (int i = 0; i < pred_horizon - 1; i++) {

		A_pow_B_cache.block(0, (i + 1)*B_.cols(), A_.rows(), B_.cols()) = (A_* A_pow_B_cache.block(0, (i)*B_.cols(), A_.rows(), B_.cols()));

	}

	du_limit(0, 0) = 70 * 1e-3 * 2;
	du_limit(0, 1) = 2;
	du_limit(1, 0) = -70 * 1e-3 * 2;
	du_limit(1, 1) = -2;
	u_limit(0, 0) = 0.6 / 2 - 0.01;
	u_limit(0, 1) = 50;
	u_limit(1, 0) = -u_limit(0, 0);
	u_limit(1, 1) = -u_limit(0, 1);



}

//When you add the const keyword to a method the this pointer will essentially become a pointer to const object, and you cannot therefore change any member data. (This is not completely true, because you can mark a member as mutable and a const method can then change it. It's mostly used for internal counters and stuff.).

Matrix<double, 2 * mpc_control_horizon, 1>& AdaGrad::Evaluate(Matrix<double, 2 * mpc_control_horizon, 1>&   x) {

	saturation_count = 0;
	if (new_desire_state) {
		new_desire_state = false;
		residuum_old = 1000.0;
		gradient_past.setZero();
	}
	x_states.block(0, 0, x0_.rows(), x0_.cols()) = x0_;
	/*std::ofstream myfile;
	myfile.open(string("jacobians")+ std::to_string(count_jacobians) + string(".csv"));
	count_jacobians++;*/
	for (int iter = 0; iter < max_iter; iter++) {
		///TODO: Ove tri linije izbacit
		//deriv_wrt_u.setZero();
		//x_states.setZero();


		u_past = 1 * u_current;



		for (int i = 0; i < pred_horizon; i++) {
			if (i < control_horizon) {
				u << x(0 * control_horizon + (i), 0), x(0 * control_horizon + (i), 0),
					x(1 * control_horizon + (i), 0), -x(1 * control_horizon + (i), 0);

			}


			x_states.block(0, i + 1, x0_.rows(), x0_.cols()) = (A_ * x_states.block(0, i, x0_.rows(), x0_.cols()) + B_ * u);
			lambdas_x.block(0, i, x0_.rows(), x0_.cols()) = -1 * x_ss_ + x_states.block(0, i, x0_.rows(), x0_.cols());



			lambdas_u.block(0, i, u_past.rows(), u_past.cols()) = u - u_past;
			lambdas_u_ref.block(0, i, u.rows(), u.cols()) = u - u_ss_;


			//derivation of u
			if (i < control_horizon) {
				//deriv_wrt_u.block(0, i, u.rows(), u.cols()) = (deriv_wrt_u.block(0, i, u.rows(), u.cols()) + (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past));
				deriv_wrt_u.block(0, i, u.rows(), u.cols()) = ((2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past));


				if (i > 0) {
					deriv_wrt_u.block(0, i - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, i - 1, u.rows(), u.cols()) - 2 * R_delta_*u);

				}
			}
			else {
				deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) + (2 * R_*u) - 2 * R_*u_ss_ + (4 * R_delta_*u) + (-2 * R_delta_*u_past));//.eval();


				deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) = (deriv_wrt_u.block(0, control_horizon - 1, u.rows(), u.cols()) - 2 * R_delta_*u);


			}
			//derivation of x
			for (int j = 0; j <= i; j++) {


				if (j < control_horizon && i < control_horizon) {

					deriv_wrt_u.block(0, j, u.rows(), u.cols()) = (deriv_wrt_u.block(0, j, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();

				}
				else {
					if (j >= control_horizon) {
						deriv_wrt_u.block(0, 4, u.rows(), u.cols()) = (deriv_wrt_u.block(0, 4, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();

					}
					else {
						deriv_wrt_u.block(0, j, u.rows(), u.cols()) = (deriv_wrt_u.block(0, j, u.rows(), u.cols()) + ((2 * Q_*x_states.block(0, i + 1, x0_.rows(), x0_.cols()) - 2 * Q_*x_ss_).transpose()*A_pow_B_cache.block(0, (i - j)*B_.cols(), A_.rows(), B_.cols())).transpose());//.eval();

					}

				}

			}

			u_past = u;
		}

		lambdas_u_ref = scale_MV_inv * lambdas_u_ref;
		lambdas_u = scale_MV_inv * lambdas_u;
		lambdas_x = scale_OV_inv * lambdas_x;




		residuum_signal = (lambdas_u_ref.cwiseProduct(R_*lambdas_u_ref)).sum() + (lambdas_u.cwiseProduct(R_delta_*lambdas_u)).sum();

		residuum_state = (lambdas_x.cwiseProduct(Q_*lambdas_x)).sum();

		residuum = residuum_signal + residuum_state;



		
		for (int iter_der = 0; iter_der < control_horizon; iter_der++) {
			//derivation df/du1
			gradients(iter_der, 0) = 2 * deriv_wrt_u(0, iter_der);
			//derivation df/du2
			gradients(iter_der + mpc_control_horizon, 0) = 2 * deriv_wrt_u(2, iter_der );
		}
		

		gradient_past = (gradient_past + gradients.cwiseProduct(gradients)).eval();


		/*if (residuum < min_residuum) {
			return x;
		}*/


		/*
		for (int print_iter = 0; print_iter < 2*mpc_control_horizon; print_iter++) {
			myfile << Jacobian( 0, print_iter) << ",";
		}
		myfile << "\n";*/


		//cout << "Jacobians" <<endl << Jacobian << endl;
		//Radi !!
		if (residuum < 0.1) {
			alfa = 1e-5;
		}

		
		for (int ada_iter = 0; ada_iter < num_heuristic_variables*mpc_control_horizon; ada_iter++) {
			change_x(ada_iter, 0) = -1 * (epsilon*gradients(ada_iter, 0)) / (sqrt(gradient_past(ada_iter, 0)) + delta);
		}

		x = (x + change_x).eval();

		/*

		if ((residuum - residuum_old) > 1e-4) {
			x = this->check_bounderies(x);
			return x;
		}
		*/

		residuum_old = residuum;

	}

	//myfile.close();
	x = this->check_bounderies(x);
	return x;
}


AdaGrad::~AdaGrad()
{
}

Matrix<double, 2 * mpc_control_horizon, 1> AdaGrad::check_bounderies(Matrix<double, 2 * mpc_control_horizon, 1>   x) {
	dummy_u = u_current;
	for (int i = 0; i < mpc_control_horizon; i++) {

		for (int j = 0; j < num_heuristic_variables; j++) {

			if (x(j * mpc_control_horizon + (i), 0) > dummy_u(2 * j, 0) + du_limit(0, j)) {
				/*cout << "dummy u " << endl << dummy_u << endl;
				cout << "du_lim" << endl << du_limit << endl;
				cout << x(j * mpc_control_horizon + (i), 0) << "   " << dummy_u(2 * j, 0) << "  " << du_limit(0, j)  << "  " << u_limit(0, j)<< endl;*/


				x(j * mpc_control_horizon + (i), 0) = dummy_u(2 * j, 0) + du_limit(0, j);




			}
			else if (x(j * mpc_control_horizon + (i), 0) < dummy_u(2 * j, 0) + du_limit(1, j)) {
				/*cout << "dummy u " << endl << dummy_u << endl;
				cout << "du_lim" << endl << du_limit << endl;
				cout << x(j * mpc_control_horizon + (i), 0) << "   " << dummy_u(2 * j, 0) << "  " << du_limit(1, j) << "  " << u_limit(1, j) << endl;*/

				x(j * mpc_control_horizon + (i), 0) = dummy_u(2 * j, 0) + du_limit(1, j);


			}

			if (x(j * mpc_control_horizon + (i), 0) > u_limit(0, j)) {
				x(j * mpc_control_horizon + (i), 0) = u_limit(0, j);
			}
			else if (x(j * mpc_control_horizon + (i), 0) < u_limit(1, j)) {
				x(j * mpc_control_horizon + (i), 0) = u_limit(1, j);
			}


		}
		//cout << "x " << endl << x << endl;
		dummy_u << x(0 * mpc_control_horizon + (i), 0), x(0 * mpc_control_horizon + (i), 0),
			x(1 * mpc_control_horizon + (i), 0), -x(1 * mpc_control_horizon + (i), 0);
	}
	return x;


}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    std::cout <<"Pokaz ---------" << nrhs << std::endl;
    /*Eigen::Matrix<double, 3,3> matricaA;
    matricaA.setZero();
    if(nrhs > 3)
    {
        
        
        
        
        
        mxAssert(mxGetClassID(prhs[2]) == mxDOUBLE_CLASS,
        "Type of the input matrix isn't double");
            mwSize     m = mxGetM (prhs[2]);
        mwSize     n = mxGetN (prhs[2]);
        mwSize    nz = mxGetNzmax (prhs[2]);
        const unsigned int a = static_cast<unsigned int>(m);
        const unsigned int c = static_cast<unsigned int>(n);    
        mexPrintf("%d\n", m);
        mexPrintf("%d\n", n);
        mexPrintf("%d\n", nz);
        std::cout <<"Pokaz ---------" << a << std::endl;
        std::cout <<"Pokaz ---------" << c << std::endl;
        matricaA = Eigen::Map<Eigen::Matrix<double, 3,3>>(mxGetPr (prhs[2]));
        std::cout << matricaA << std::endl;
        //std::cout <<"Pokaz ---------" << proba  << std::endl;
    }*/
    
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
        
    // New
    if (!strcmp("new", cmd)) {
        Matrix<double, num_state_variables, num_state_variables> A;
        Matrix<double, num_state_variables, num_manipulated_variables>  B;
        Matrix<double, num_state_variables, num_state_variables> Bd;
        Matrix<double, num_state_variables, num_state_variables> Q;
        Matrix<double, num_state_variables, num_state_variables> Q_final;
        Matrix<double, num_manipulated_variables, num_manipulated_variables> R;
        Matrix<double, num_manipulated_variables, num_manipulated_variables> R_delta;
        Matrix<double, num_state_variables, 1> disturbance;
        int num_params;
        int pred_horizon;
        int control_horizon;
        Matrix<double, num_manipulated_variables, num_manipulated_variables> scale_MV;
        Matrix<double, num_state_variables, num_state_variables> scale_OV;
        
        A.setZero();
        B.setZero();
        Bd.setZero();
        Q.setZero();
        Q_final.setZero();
        R.setZero();
        R_delta.setZero();
        disturbance.setZero();
        scale_MV.setZero();
        scale_OV.setZero();
        
        A = Eigen::Map<Matrix<double, num_state_variables, num_state_variables>>(mxGetPr (prhs[1]));
        B = Eigen::Map<Matrix<double, num_state_variables, num_manipulated_variables>>(mxGetPr (prhs[2]));
        Bd = Eigen::Map<Matrix<double, num_state_variables, num_state_variables>>(mxGetPr (prhs[3]));
        Q = Eigen::Map<Matrix<double, num_state_variables, num_state_variables>>(mxGetPr (prhs[4]));
        Q_final = Eigen::Map<Matrix<double, num_state_variables, num_state_variables>>(mxGetPr (prhs[5]));
        R = Eigen::Map<Matrix<double, num_manipulated_variables, num_manipulated_variables>>(mxGetPr (prhs[6]));
        R_delta = Eigen::Map<Matrix<double, num_manipulated_variables, num_manipulated_variables>>(mxGetPr (prhs[7]));
        disturbance = Eigen::Map<Matrix<double, num_state_variables, 1>>(mxGetPr (prhs[8]));
        num_params = (int)mxGetScalar(prhs[9]);
        pred_horizon = (int)mxGetScalar(prhs[10]);
        control_horizon = (int)mxGetScalar(prhs[11]);
        scale_MV = Eigen::Map<Matrix<double, num_manipulated_variables, num_manipulated_variables>>(mxGetPr (prhs[12]));
        scale_OV = Eigen::Map<Matrix<double, num_state_variables, num_state_variables>>(mxGetPr (prhs[13]));
        
        //std::cout << num_params << std::endl;
        
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        // Return a handle to a new C++ instance
        
        plhs[0] = convertPtr2Mat<AdaGrad>(new AdaGrad(A, B, Bd, Q, Q_final, R, R_delta, disturbance,num_params, pred_horizon,control_horizon,scale_MV, scale_OV  ));
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<AdaGrad>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    AdaGrad* AdaGrad_instance = convertMat2Ptr<AdaGrad>(prhs[1]);
    
    // Call the various class methods
    // Train    
    if (!strcmp("Evaluate", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Train: Unexpected arguments.");
        
        // Call the method
        
        
        Matrix<double, 2 * mpc_control_horizon, 1> x_ =  Eigen::Map<Matrix<double, 2 * mpc_control_horizon, 1>>(mxGetPr (prhs[2]));
        //std::cout << "AJMOOOO " << nrhs << std::endl;
        //x_.setZero();
        Matrix<double, 2 * mpc_control_horizon, 1> resultEigen = AdaGrad_instance->Evaluate(x_);
        //Matrix<double, 2 * mpc_control_horizon, 1> resultEigen ;
        //resultEigen.setZero();
        //Matrix<double, 2 * mpc_control_horizon, 1> resultEigen;
        //resultEigen.setZero();
        mwSize rows = resultEigen.rows();
        mwSize cols = resultEigen.cols();
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL); // Create MATLAB array of same size
        Eigen::Map<Matrix<double, 2 * mpc_control_horizon, 1>> map(mxGetPr(plhs[0]), rows, cols); // Map the array
        map = resultEigen; // Copy
        
        return;
        
    }
    
      
    if (!strcmp("set_current_state", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Train: Unexpected arguments.");
        
        // Call the method
        
        
        Matrix<double, num_state_variables, 1> x_ =  Eigen::Map<Matrix<double,num_state_variables, 1>>(mxGetPr (prhs[2]));
        AdaGrad_instance->set_x0_(x_);
        
        
        return;
        
    }
    
    if (!strcmp("set_desire_state", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Train: Unexpected arguments.");
        
        // Call the method
        
        
        Matrix<double, num_state_variables, 1> x_ =  Eigen::Map<Matrix<double,num_state_variables, 1>>(mxGetPr (prhs[2]));
        AdaGrad_instance->set_x_ss(x_);
        
        
        return;
        
    }
   
    
    if (!strcmp("set_desire_stacionary_input", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Train: Unexpected arguments.");
        
        // Call the method
        
        
        Matrix<double, num_manipulated_variables, 1> x_ =  Eigen::Map<Matrix<double,num_manipulated_variables, 1>>(mxGetPr (prhs[2]));
        AdaGrad_instance->set_u_ss(x_);
        
        
        return;
        
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
