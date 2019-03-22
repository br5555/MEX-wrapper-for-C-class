#ifndef ADA_GRAD_MPC_SOLVER
#define ADA_GRAD_MPC_SOLVER

#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <string> 


using namespace Eigen;
using namespace std;

/// <summary>
/// defines number of state variables in linearized model
/// </summary>
constexpr int num_state_variables = 8;
/// <summary>
/// defines how many input variables we have 
/// </summary>
constexpr int num_manipulated_variables = 4;
/// <summary>
/// how many DOF we really have ( using symmetry u1 = u2 and u3 = -u4 )
/// </summary>
constexpr int num_heuristic_variables = 2; //using symmetry u1 = u2 and u3 = -u4 
/// <summary>
/// MPC parameter which tells how many prediction horizon steps are optimized
/// </summary>
constexpr int mpc_control_horizon = 5;
/// <summary>
/// MPC parameter prediction horizon which tells how many steps in the future are unroll
/// </summary>
constexpr int prediction_horizon = 14;
/// <summary>
/// parameter defines number of cost functions or number of criteria
/// </summary>
constexpr int num_cost_functions = 1;
class AdaGrad
{
public:
	/// <summary>
	/// 
	/// </summary>
	/// <param name="A"></param>
	/// <param name="B"></param>
	/// <param name="Bd"></param>
	/// <param name="Q"></param>
	/// <param name="Q_final"></param>
	/// <param name="R"></param>
	/// <param name="R_delta"></param>
	/// <param name="disturbance"></param>
	/// <param name="num_params"></param>
	/// <param name="pred_horizon"></param>
	/// <param name="control_horizon"></param>
	/// <param name="scale_MV"></param>
	/// <param name="scale_OV"></param>
	explicit AdaGrad(Matrix<double, num_state_variables, num_state_variables> A, Matrix<double, num_state_variables, num_manipulated_variables>  B,Matrix<double, num_state_variables, num_state_variables>Bd,Matrix<double, num_state_variables, num_state_variables> Q,Matrix<double, num_state_variables, num_state_variables> Q_final, Matrix<double, num_manipulated_variables, num_manipulated_variables> R,Matrix<double, num_manipulated_variables, num_manipulated_variables> R_delta,Matrix<double, num_state_variables, 1> disturbance, int num_params, int pred_horizon, int control_horizon,Matrix<double, num_manipulated_variables, num_manipulated_variables> scale_MV, Matrix<double, num_state_variables, num_state_variables> scale_OV);
	
	/// /// <summary>
	/// 
	/// </summary>
	~AdaGrad();
	/// <summary>
	/// 
	/// </summary>
	/// <param name="x"></param>
	/// <returns></returns>
	Matrix<double, 2 * mpc_control_horizon, 1>&  Evaluate(Matrix<double, 2 * mpc_control_horizon, 1>&   x);
	/// <summary>
	/// 
	/// </summary>
	/// <returns></returns>
	int dim_X(void) const { return 12; }
	void set_u_past(Matrix<double, num_manipulated_variables, 1> u_past_) { this->u_past = u_past_; }
	void set_u_current(Matrix<double, num_manipulated_variables, 1> u_current_) { this->u_current = u_current_; }
	void set_u_ss(Matrix<double, num_manipulated_variables, 1> u_ss) { this->u_ss_ = u_ss; }
	void set_x_ss(Matrix<double, num_state_variables, 1>  x_ss) { this->x_ss_ = x_ss; new_desire_state = true; }
	void set_x0_(Matrix<double, num_state_variables, 1>  x0) { this->x0_ = x0; }
	void set_A(Matrix<double, num_state_variables, num_state_variables> A) { this->A_ = A; }
	void set_B(Matrix<double, num_state_variables, num_manipulated_variables> B) { this->B_ = B; }
	void set_Bd(Matrix<double, num_state_variables, num_state_variables> Bd) { this->Bd_ = Bd; }
	void set_Q(Matrix<double, num_state_variables, num_state_variables>  Q) { this->Q_ = Q; }
	void set_Q_final(Matrix<double, num_state_variables, num_state_variables>  Q_final) { this->Q_final_ = Q_final; }
	void set_R(Matrix<double, num_manipulated_variables, num_manipulated_variables>  R) { this->R_ = R; }
	void set_R_delta(Matrix<double, num_manipulated_variables, num_manipulated_variables>  R_delta) { this->R_delta_ = R_delta; }
	void set_insecure(Matrix<double, num_state_variables, 1> insecure) { this->insecure_ = insecure; }
	void set_disturbance(Matrix<double, num_state_variables, 1> disturbance) {
		this->disturbance_ = disturbance;
		this->insecure_ = this->Bd_*this->disturbance_;
	}
	void set_num_params(int num_params) { this->num_params_ = num_params; }
	double get_residuum() { return residuum; }
	double get_residuum_signal() { return residuum_signal; }
	double get_residuum_state() { return residuum_state; }
	MatrixXd get_x_states() { return x_states; }
	/*If you define a structure having members of fixed - size vectorizable Eigen types, you must overload its "operator new" so that it generates 16 - bytes - aligned pointers.*/
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	Matrix<double, 2 * mpc_control_horizon, 1>  check_bounderies(Matrix<double, 2 * mpc_control_horizon, 1>   x);
	bool new_desire_state;
	int  num_params_, pred_horizon, control_horizon, max_iter, saturation_count, count_jacobians;
	double residuum, residuum_signal, residuum_state, epsilon, delta, min_residuum, alfa, residuum_old;
	Matrix<double, num_state_variables, num_state_variables> A_, Bd_;
	Matrix<double, num_state_variables, num_manipulated_variables> B_;
	Matrix<double, num_manipulated_variables, num_manipulated_variables> R_, R_delta_, scale_MV_inv;
	Matrix<double, num_state_variables, num_state_variables> Q_, Q_final_, scale_OV_inv;
	Matrix<double, num_state_variables, 1> x_ss_, disturbance_, insecure_, x0_;
	Matrix<double, num_manipulated_variables, 1> u_ss_, u_prev_, u, u_past, u_current, dummy_u;
	Matrix<double, mpc_control_horizon*num_heuristic_variables, num_cost_functions> gradients, gradient_past;
	Matrix<double, num_state_variables, prediction_horizon * num_manipulated_variables> A_pow_B_cache;
	Matrix<double, num_state_variables, prediction_horizon>  lambdas_x;
	Matrix<double, num_manipulated_variables, prediction_horizon>  lambdas_u, lambdas_u_ref;
	Matrix<double, num_state_variables, prediction_horizon + 1>  x_states; // + 1 x0 is added 
	Matrix<double, num_manipulated_variables, mpc_control_horizon> deriv_wrt_u;
	Matrix<double, 2 * mpc_control_horizon, 1>  change_x;
	Matrix<double, 2, num_heuristic_variables>  du_limit, u_limit; //2 for upper and lower





};
#endif 