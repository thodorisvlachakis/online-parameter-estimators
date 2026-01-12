function dXdt = gradient_descent_online_estimator_MSD(t, X, m, b, k, u, Gamma, a_m)
% This function considers the mass-spring-damper system with an external force,
% which is described by the differential equation: m*x''(t) + b*x'(t) + k*x(t) = u(t),
% where x(t) is the is the displacement from the equilibrium position, m>0
% the mass, k>0 the lattice constant, b>0 a constant damping coefficient and
% u(t) the external force. The function defines the state equations of this
% linear system considering the parameters m, k and b of the system as known
% and then assumes these parameters as unknown and defines a design of real-time
% estimators for them, exploiting the measurable states of the real system
% and the measurable input u.
% Main Idea: Thε function gradient_descent_online_estimator_MSD implements
% the design of a real-time estimator based on the gradient descent method
% as an online estimation method (i.e adaptive estimation), for the estimation
% of the unknown parameters k, m and b of the mass-spring-damper system. The
% real system is considered as available through its state variables and the
% input u, which are measurable, and so the function designs a model for 
% estimating the real system, through the estimation of its unknown parameters
% with gradient descent (online) method. In that way, it also calculates
% estimations for the state variables of the real system.
% The gradient descent method is applied to systems that can be expressed in
% linearly parameterized form. For the mass-spring-damper system, the second
% equation of state is as follows: x2' = -theta1*x1 - theta2*x2 + theta3*u,
% where theta1 = k/m, theta2=b/k and theta3 = 1/m. Therefore, the system can
% be expressed in linearly parameterized form as follows: x2 = (theta_star)^T * phi,
% where phi is a vector whose components include only the measurable states
% x1, x2 of the system and the measurable input u and theta_new is a parameter
% vector whose components include the (composite) parameters theta1, theta2,
% theta3. Finally, The function gradient_descent_online_estimator_MSD implements
% an online estimator for the parameter vector theta_star.
%
% The function takes as input the following:
% - t: the simulation time.
% - X: an extended vector of length 8 whose first 2 elements represent the 
%      state vector of the system, the next 1 element represent the estimation
%      of first state variable of the real system, the next 3 elements 
%      represent the components of the vector phi -that includes the
%      information from the measurements of state vector and input u- of
%      the linearly parameterized expression x2 = (theta_star)^T * phi that
%      describes the system and the last 3 elements represents the real-time 
%      estimator for the unknown parameters theta_star1 = k/m, theta_star2=b/k - a_m
%      and theta_star3 = 1/m.
% - m: a positive number which represents the mass.
% - k: a positive number which represents the lattice constant.
% - b: a positive number which represents the constant damping coefficient.
% - u: a function of time that represents the external force (input) of the
%      system.
% - Gamma: a diagonal matrix whose non-zero element of each row represents
%          the gain (learning-rate) of the corresponing real-time estimator
%          for each composite parameter theta_star1= k/m, theta_star2= b/k - a_m,
%          theta_star3 = 1/m.
% - a_m: a positive number that serves the definition of a stable first-order
%        filter that is used to bring the system into a linearly parameterized
%        form. The a_m appears in the dynamics of the vector phi with the
%        "filtered" measurable signals.
    
    % Decoding of the extended state vector
    x = X(1:2);
    %x1_hat = X(3);
    phi = X(4:6);
    theta_star_hat = X(7:9); % theta_star1_hat = X(8); theta_star2_hat = X(9); theta_star3_hat = X(10);
    
    % Real mass-spring-damper system with an external force
    A = [0 1;
        -k/m -b/m];
    
    B = [0;
        1/m];
    
    dxdt = systemEquationsOfState(t, x, A, B, u);
    
    % Define the estimation of the second state variable exploiting the 
    % linearly parameterized form and then compute the modeling error of
    % the second variable state. Also design the dynamics for the estimation
    % of the first state variable of the real system, using a model that
    % considers the stucture of the real system.
    x2_hat = theta_star_hat' * phi ;
    
    error_x2 = x(2) - x2_hat;
    
    dx1dt_hat = x2_hat;
    
    % Real-time estimator
    if( isa(u, 'function_handle') )
        % if u is a function of time compute its value for the specific time t.
        u_val = u(t);
    else
        u_val = u;
    end
    
    dtheta_stardt_hat = Gamma * error_x2 * phi;
    dphidt = -a_m * phi + [-x(1); -x(2); u_val];

    % Design the model estimating the real system
    
    % Dynamics of the extended vector of states and estimations, Χ
    dXdt = [dxdt; dx1dt_hat; dphidt; dtheta_stardt_hat];
end