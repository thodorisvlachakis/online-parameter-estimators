function dXdt = composite_structure_model_adaptive_estimation_RollAngleSys(t, X, a1, a2, a3, b, u, Gamma, theta_m, d)
% This function considers the roll angle system with control input, which is
% is described by the following non-linear differential equation:
% r''(t) = - a1 * r'(t) - a2 * sin(r(t)) + a3 * (r'(t))^2 * sin(2*r(t)) + b * u(t) + d(t),
% where r(t) [rad] is the roll angle, a1,a2,a3 > 0 and b > 0 are constant
% parameters of the dynamic system, u(t) is the system control input and
% d(t) is external disturbances.
% The function defines the state equations of the above system considering 
% the parameters a1, a2, a3 and b of this system as known and then assumes
% these parameters as unknown and defines a design of real-time estimators 
% for them, exploiting the measurable states of the real system and the 
% measurable input u and considering that the nonlinear functions that appear
% in the above linear differential equation of the roll angle system are known.
% Main Idea: Thε function composite_model_structure_adaptive_estimation_RollAngleSys
% implements the design of a real-time estimator of a composite structure 
% (topology) based on the Lyapunov method for adaptive estimation, for the
% estimation of the unknown parameters a1, a2, a3 and b of the roll angle 
% system. The real system is considered as available through its state variables
% and the input u, which are measurable, and so the function designs a model 
% for estimating the real system, using the composite structure based on adaptive
% parameter estimation, through the Lyapunov method.
%
% The function takes as input the following:
% - t: the simulation time.
% - X: an extended vector of length 8 whose first 2 elements represent the 
%      state vector of the system, the next 2 elements represent the
%      estimations of the real system's state variables (respectively) and
%      the last 4 elements represents the real-time estimator for the
%      unknown parameters a1, a2, a3 and b.
% - a1: a positive number which represents the parameter a1
% - a2: a positive number which represents the parameter a2
% - a3: a positive number which represents the parameter a3
% - b: a positive number which represents the parameter b
% - u: a function of time that represents the control input of the system
%      or a state feedback controller, which is actually a function of both
%      time t and state vector x.
% - Gamma: a diagonal matrix whose non-zero element of each row represents
%          the gain (learning-rate) of the corresponing real-time estimator
%          for each parameter a1, a2, a3 and b.
% - theta_m: a vector of length 2 whose each element represents the corresponding
%            model following gains for each state variable of the real
%            system (gains of the corresponding modeling gains in the expression
%            of the composite structure model).
% - d: a function of time or a vector representing the possible external
%      disturbances
    
    % Declare if there are external disturbances or not.
    if nargin<10
        d = 0;
    end

    % Decoding of the extended state vector
    x = X(1:2);
    x_hat = X(3:4);
    theta_hat = X(5:8); % a1_hat = X(5); a2_hat = X(6); a3_hat = X(7); b_hat = X(8)
    
    % Roll angle dynamic system
    dxdt = rollAngleDynamicSystem(t, x, a1, a2, a3, b, u, d);

    % Composite Model Structure
    if( isa(u, 'function_handle') )
        if nargin(u) == 1
            % u is a function of time, so compute its value for the specific time t.
            u_val = u(t);
        elseif nargin(u) == 2
            % u is state feedback control input
            if length(u(t,x)) > 1
                [u_val, ~ , ~] = u(t,x);
            else
                u_val = u(t,x);
            end
        end
    else
        u_val = u;
    end

    phi = [-x(2); -sin(x(1)); (x(2)^2 * sin(2*x(1))); u_val];
    
    dxdt_hat = zeros(2,1);

    dxdt_hat(1) = x(2) + theta_m(1) * (x(1) - x_hat(1));
    dxdt_hat(2) = theta_hat' * phi + theta_m(2) * (x(2) - x_hat(2));
    
    % Define modeling error of the second state variable
    error_x2 = x(2) - x_hat(2);
    
    % Real-time estimator
    dthetadt_hat = Gamma * error_x2 * phi;
    
    % Dynamics of the extended vector of states and estimations, Χ
    dXdt = [dxdt; dxdt_hat; dthetadt_hat];
end