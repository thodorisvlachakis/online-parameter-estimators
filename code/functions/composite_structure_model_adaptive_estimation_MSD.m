function dXdt = composite_structure_model_adaptive_estimation_MSD(t, X, m, b, k, u, Gamma, theta_m, n)
% This function considers the mass-spring-damper system with an external force,
% which is described by the differential equation: m*x''(t) + b*x'(t) + k*x(t) = u(t),
% where x(t) is the is the displacement from the equilibrium position, m>0
% the mass, k>0 the lattice constant, b>0 a constant damping coefficient and
% u(t) the external force. The function defines the state equations of this
% linear system considering the parameters m, k and b of the system as known
% and then assumes these parameters as unknown and defines a design of real-time
% estimators for them, exploiting the measurable states of the real system
% and the measurable input u.
% Main Idea: Thε function composite_model_structure_adaptive_estimation_MSD
% implements the design of a real-time estimator of a composite structure 
% (topology) based on the Lyapunov method for adaptive estimation, for the
% estimation of the unknown parameters k, m and b of the mass-spring-damper 
% system. The real system is considered as available through its state variables
% and the input u, which are measurable, and so the function designs a model 
% for estimating the real system, using the composite structure based on adaptive
% parameter estimation, through the Lyapunov method.
%
% The function takes as input the following:
% - t: the simulation time.
% - X: an extended vector of length 7 whose first 2 elements represent the 
%      state vector of the system, the next 2 elements represent the
%      estimations of the real system's state variables (respectively) and
%      the last 3 elements represents the real-time estimator for the
%      unknown parameters theta1 = k/m, theta2 = b/k and theta3 = 1/m.
% - m: a positive number which represents the mass.
% - k: a positive number which represents the lattice constant.
% - b: a positive number which represents the constant damping coefficient.
% - u: a function of time that represents the external force (input) of the
%      system.
% - Gamma: a diagonal matrix whose non-zero element of each row represents
%          the gain (learning-rate) of the corresponing real-time estimator
%          for each composite parameter theta1= k/m, theta2= b/k,
%          theta3 = 1/m.
% - theta_m: a vector of length 2 whose each element represents the corresponding
%            model following gains for each state variable of the real
%            system.
% - n: a function of time or a constant representing the possible noise in
%      the measurement of the output (which is the first state variable
%      x1(t) of the real system) of the real system.
    
    % Declare if there is noise in the measurement of the output x1(t)=x(t)
    % or not.
    if nargin<9
        n = 0;
    end

    if( isa(n, 'function_handle') )
        % if n is a function of time compute its value for the specific time t.
        n_val = n(t);
    else
        n_val = n;
    end

    % Decoding of the extended state vector
    x = X(1:2);
    x_hat = X(3:4);
    theta_hat = X(5:7); % theta1_hat = X(5); theta2_hat = X(6); theta3_hat = X(7);
    
    % Real mass-spring-damper system with an external force
    A = [0 1;
        -k/m -b/m];
    
    B = [0;
        1/m];
    
    dxdt = systemEquationsOfState(t, x, A, B, u);
    
    % State vector measurements
    x_measured = x + [n_val; 0] ;

    % Composite Model Structure
    if( isa(u, 'function_handle') )
        % if u is a function of time compute its value for the specific time t.
        u_val = u(t);
    else
        u_val = u;
    end
    phi = [-x_measured(1); -x_measured(2); u_val];
    
    dxdt_hat = zeros(2,1);

    dxdt_hat(1) = x_measured(2) + theta_m(1) * (x_measured(1) - x_hat(1));
    dxdt_hat(2) = theta_hat' * phi + theta_m(2) * (x_measured(2) - x_hat(2));
    
    % Define modeling error of the second state variable
    error_x2 = x_measured(2) - x_hat(2);
    
    % Real-time estimator
    dthetadt_hat = Gamma * error_x2 * phi;
    
    % Dynamics of the extended vector of states and estimations, Χ
    dXdt = [dxdt; dxdt_hat; dthetadt_hat];
end