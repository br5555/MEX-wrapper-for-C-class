clc
clear all
close all
%function run_example()
% Check if the mex exists
dir = fileparts(mfilename('fullpath'));
if ~isequal(fileparts(which('example_mex')), dir)
    % Compile the mex
    cwd = cd(dir);
    cleanup_obj = onCleanup(@() cd(cwd));
    fprintf('Compiling example_mex\n');
    mex example_mex.cpp
end



combined_control_mpc_use_ = 1;  

kStateSize = 8;
kInputSize = 4;
kMeasurementSize = 1;                       
kDisturbanceSize = kStateSize;              

kControlHorizonSteps = 5;
kPredictionHorizonSteps = 14;




	OV_scale = zeros(8,8);
	OV_scale(1, 1) = 0.58;
	OV_scale(2, 2) = 4.0;
	OV_scale(3, 3) = 0.58;
	OV_scale(4, 4) = 4.0;
	OV_scale(5, 5) = 800;
	OV_scale(6, 6) = 800;
	OV_scale(7, 7) = 0.5236;
	OV_scale(8, 8) = 0.5236;

    
    MV_scale = zeros(4,4);
	MV_scale(1, 1) = 0.58;
	MV_scale(2, 2) = 0.58;
	MV_scale(3, 3) = 100;
	MV_scale(4, 4) = 100;







    R = zeros(4,4);
	R(1, 1) = 0.67670;
	R(2, 2) = 0.67670;
	R(3, 3) = 0.13534000;
	R(4, 4) = 0.13534000;

    R_delta = zeros(4,4);
	R_delta(1, 1) = 0.738879858135067;
	R_delta(2, 2) = 0.738879858135067;
	R_delta(3, 3) = 0.007388798581351;
	R_delta(4, 4) = 0.007388798581351;


    Q = zeros(8,8);
	Q(1, 1) = 0.135340000000000;
	Q(2, 2) = 0.002706800000000;
	Q(3, 3) = 0.1353400;
	Q(4, 4) = 0.002706800;
	Q(5, 5) = 0.002706800;
	Q(6, 6) = 0.002706800;
	Q(7, 7) = 10.7068000;
	Q(8, 8) = 9.676700;

	Q_final = Q;





	

	model_Bd_ = [0.0368517, 0.000550615, 0, 0, 0, 0, 0, 0;
		-0.215525, 0.0225941, 0, 0, 0, 0, 0, 0;
		0, 0, 0.0368517, 0.000550615, 0, 0, 0, 0;
		0, 0, -0.215525, 0.0225941, 0, 0, 0, 0;
		0, 0, 0, 0, 0.0369936, 0, 0, 0;
		0, 0, 0, 0, 0, 0.0369936, 0, 0;
		4.15876e-05, 1.98561e-06, 4.15876e-05, 1.98561e-06, 4.2151e-07, -4.2151e-07, 0.04, 0.000792;
		0.00294716, 0.000146518, 0.00294716, 0.000146518, 3.13605e-05, -3.13605e-05, 0, 0.04];

	estimated_disturbances_ = [ 0.0744648;
		-0.00684855;
		-0.468717;
		0.553529;
		20;
		20;
		-0.139467;
		-0.00190078];
	model_A_70_ms =[ 0.4271, 0.0223, 0, 0, 0, 0, 0, 0;
		-8.7243, -0.1501, 0, 0, 0, 0, 0, 0;
		0, 0, 0.4271, 0.0223, 0, 0, 0, 0;
		0, 0, -8.7243, -0.1501, 0, 0, 0, 0;
		0, 0, 0, 0, 0.7300, 0, 0, 0;
		0, 0, 0, 0, 0, 0.7300, 0, 0;
		0.0090, 0.0005, 0.0090, 0.0005, 0.0001, -0.0001, 1.0000, 0.0787;
		0.1699, 0.0113, 0.1699, 0.0113, 0.0028, -0.0028, 0, 1.0000];

	model_B_70_ms =[ 0.5729, 0, 0, 0;
		8.7243, 0, 0, 0;
		0, 0.5729, 0, 0;
		0, 8.7243, 0, 0;
		0, 0, 0.2700, 0;
		0, 0, 0, 0.2700;
		-0.0037, -0.0037, 0.0000, -0.0000;
		-0.0347, -0.0347, 0.0005, -0.0005];





% Use the example interface
% This is the interface written specifically for the example class
fprintf('Using the example interface\n');
obj = example_interface(model_A_70_ms, model_B_70_ms, model_Bd_, Q, Q_final, R, R_delta, estimated_disturbances_, kStateSize, 14, kControlHorizonSteps, MV_scale, OV_scale);

u_input = zeros(10,1);

x_ss = zeros(8,1);
x_ss(7,1) = 15*(pi/180);

u_ss = zeros(4,1);

x = zeros(8,1);

display('before calling')

set_desire_state(obj, x_ss);
input_real = zeros(4,1);
set_desire_stacionary_input(obj, u_ss);

for i = 1 : 100
set_current_state(obj, x);
u_input = Evaluate(obj, u_input);
input_real(1,1) = u_input(1,1);
input_real(2,1) = u_input(1,1);
input_real(3,1) = u_input(kControlHorizonSteps,1);
input_real(4,1) = -u_input(kControlHorizonSteps,1);
x = model_A_70_ms*x + model_B_70_ms * input_real;

end
display('after calling')
clear obj % Clear calls the delete method

% Use the standard interface
% This interface can be used for any mex interface function using the
% pattern:
%   Construction -    obj = mexfun('new',         ...)
%   Destruction -           mexfun('delete', obj)
%   Other methods - [...] = mexfun('method', obj, ...)
% The standard interface avoids the need to write a specific interface
% class for each mex file.



fprintf('Using the standard interface\n');
obj = mex_interface(str2fun([dir '/example_mex']), model_A_70_ms, model_B_70_ms, model_Bd_, Q, Q_final, R, R_delta, estimated_disturbances_, kStateSize, 14, kControlHorizonSteps, MV_scale, OV_scale); % str2fun allows us to use the full path, so the mex need not be on our path
for i = 1 : 100
obj.set_current_state(x);
u_input = obj.Evaluate(u_input);
input_real(1,1) = u_input(1,1);
input_real(2,1) = u_input(1,1);
input_real(3,1) = u_input(kControlHorizonSteps,1);
input_real(4,1) = -u_input(kControlHorizonSteps,1);
x = model_A_70_ms*x + model_B_70_ms * input_real;

end
clear obj % Clear calls the delete method
