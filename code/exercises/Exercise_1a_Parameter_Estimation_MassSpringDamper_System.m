% Task a) below.

% ------------------------ Exercise 1 ------------------------
% In this exercise we consider the mass-spring-damper system with an external
% force, which is described by the differential equation: m*x''(t) + b*x'(t) + k*x(t) = u(t),
% where x(t) is the is the displacement from the equilibrium position, m>0
% the mass, k>0 the lattice constant, b>0 a constant damping coefficient and
% u(t) the external force. In this exercise we will consider such a real system,
% which is available for measurement (meaning that the vector of the system's
% state variables, as well as the external force u(t) that acts as an input
% to the system, are measurable signals) and we will implement real-time 
% estimators, for estimating the unknown parameters m, k and b. 
% In the following tasks, we will perform a simulation of the mass-spring-damper
% system using an ode solver in order to simulate the response of this system.
% The simulation of the system will make available the state vector, which
% we assume is available for measurement in real-time parameter estimation
% methods. Obviously, simulating the system requires knowledge of it, i.e. 
% the actual values ​​of the parameters m, k and b. For our experiments, we 
% will define the actual values ​​of m, k, b and design real-time estimators 
% of these parameters, using first the gradient method (task a) and then the
% Lyapunov method (task b) - using both the parallel topology (b1) and the 
% composite topology (b2). In each of the following estimator designs, we 
% initially define the estimation model of the real system and choose the
% initial conditions for this model, as well as an initial estimations for
% the parameters (which acts as the initial condition for the real-time estimator).

% * In each of the following tasks we design real-time estimators for a
% parameter vector theta = [k/m b/m 1/m]^T and then we can extract the
% real-tme estimators for the exact parameters m, k and b, just considering
% the elements of the vector theta.

% ** In the following tasks, the real-time estimator is designed with the 
% corresponding method each time and a simulation of the real system and the
% model, designed for the estimation of the real system, is performed for a
% simulation time period of 20 sec with an appropriate integration step for
% accurate results.

addpath(fullfile(pwd, '..', 'functions'));

% Declare the actual values of the parameters m, k and b
m_real = 1.315;
k_real = 0.725;
b_real = 0.225;

% Declare the initial conditions of the system
x0 = [0; 0];

% Define the simulation time to be 20 sec and declare a strict condtion
% for the integration step Δt in order to achieve more accurate results.
tspan = [0 20];
integration_step_lim = 1e-3;

options = odeset('MaxStep',integration_step_lim);

figno = 0;

% a) In this task we will design a real-time estimator of the unknown parameters
% m, k and b of the real system, using the gradient descent (on-line)
% method. We will perform the experiment twice. One with a constant input u
% and the other with a sinusoidal.
% When using the gradient descent method, it is necessary to express the 
% system in linearly parameterized form. To define the vector (of time t) phi,
% whose components include the measurable signals, we need to define the 
% parameters of a stable filter with which we essentially filter the measurable
% signals. The filtered outputs are the components of phi. The states x(t),
% x'(t), and the input u of the system are considered measurable signals. 
% The filter we will need is a first-order one.
% Note that the model we design The model determines the estimation for the
% second state variable of the system (x2) from the linearly parameterized
% expression x2_hat = theta_star_hat^T * phi, which uses the parameter 
% estimators, while the first state variable (x1) is determined by the dynamic
% (x1_hat)'(t) = x2_hat(t) (and an initial condition for x1_hat) that we 
% define in the model based on knowledge of the system structure.

% Define the parameter a_m of the filter.
a_m = 2;

% Declare the initial conditions for the model of the real system (only the
% initial condition for the first state variable is needed), the 
% measurements in phi and the real-time estimator
x0_hat = 0;
phi0 = [0; 0; 0];
theta_star0_hat = [2; 0.2; 0.5];

% Declare the initial conditions for the extended state vector X which will
% be used to perform the simulation of the real system and compute the 
% and the real-time estimators for the parameters.
X0 = [x0; x0_hat; phi0; theta_star0_hat];

% a1) Case 1: Constant Input u.
% Define the input of the system
u = 2.5 ;

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta_star = [k/m (b/m - a_m) 1/m]^T
Gamma_grad_desc_1 = diag([1, 1.1, 0.095]);

[t_grad_desc_1, X_grad_desc_1] = ode45(@(t,X) gradient_descent_online_estimator_MSD(t,X,m_real,b_real,k_real,u,Gamma_grad_desc_1,a_m), tspan, X0, options);

x_grad_desc_1 = X_grad_desc_1(:, 1:2);
phi_grad_desc_1 = X_grad_desc_1(:, 4:6);
theta_star_hat_grad_desc_1 = X_grad_desc_1(:, 7:9);
%x2_hat_grad_desc_1 = sum(phi_grad_desc_1.*theta_star_hat_grad_desc_1 , 2);
x_hat_grad_desc_1 = cat(2, X_grad_desc_1(:, 3), sum(phi_grad_desc_1.*theta_star_hat_grad_desc_1 , 2) );

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_star_hat is the estimation of 
% theta_star = [k/m (b/m - a_m) 1/m]^T and solving the equations resulting
% from the corresponding components.
m_hat_grad_desc_1 = 1./theta_star_hat_grad_desc_1(:,3);
k_hat_grad_desc_1 = m_hat_grad_desc_1.*theta_star_hat_grad_desc_1(:,1);
b_hat_grad_desc_1 = m_hat_grad_desc_1.*(theta_star_hat_grad_desc_1(:,2) + a_m);

final_estimations_grad_desc_1 = [m_hat_grad_desc_1(end), k_hat_grad_desc_1(end), b_hat_grad_desc_1(end)];

% Create the graphical representations of the output x(t) (first component
% of the state vector x of the system) of the real system, the estimate x_hat(t)
% from the designed model, as well as the difference error_x(t) = x(t) - x_hat(t)
% which represents the modeling error of the output. Additionally, create
% the graphical representations of the real-time estimators m_hat(t), k_hat(t)
% and b_hat(t) for the time the simulation was run.

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t_grad_desc_1, x_grad_desc_1(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_grad_desc_1, x_hat_grad_desc_1(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_grad_desc_1, x_grad_desc_1(:,1)-x_hat_grad_desc_1(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Gradient Descent - Constant Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_grad_desc_1, m_hat_grad_desc_1, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_1, m_real*ones(size(t_grad_desc_1)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Gradient Descent - Constant Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_grad_desc_1, k_hat_grad_desc_1, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_1, k_real*ones(size(t_grad_desc_1)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Gradient Descent - Constant Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_grad_desc_1, b_hat_grad_desc_1, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_1, b_real*ones(size(t_grad_desc_1)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Gradient Descent - Constant Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on


% a2) Case 2: Sinusoidal Input u.
% Define the input of the system
u = @(t) 2.5*sin(t) ;

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta_star = [k/m (b/m - a_m) 1/m]^T
Gamma_grad_desc_2 = diag([0.5, 8, 1.5]);

[t_grad_desc_2, X_grad_desc_2] = ode45(@(t,X) gradient_descent_online_estimator_MSD(t,X,m_real,b_real,k_real,u,Gamma_grad_desc_2,a_m), tspan, X0, options);

x_grad_desc_2 = X_grad_desc_2(:, 1:2);
phi_grad_desc_2 = X_grad_desc_2(:, 4:6);
theta_star_hat_grad_desc_2 = X_grad_desc_2(:, 7:9);
%x2_hat_grad_desc_2 = sum(phi_grad_desc_2.*theta_star_hat_grad_desc_2 , 2);
x_hat_grad_desc_2 = cat(2, X_grad_desc_2(:, 3), sum(phi_grad_desc_2.*theta_star_hat_grad_desc_2 , 2) );

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_star_hat is the estimation of 
% theta_star = [k/m (b/m - a_m) 1/m]^T and solving the equations resulting
% from the corresponding components.
m_hat_grad_desc_2 = 1./theta_star_hat_grad_desc_2(:,3);
k_hat_grad_desc_2 = m_hat_grad_desc_2.*theta_star_hat_grad_desc_2(:,1);
b_hat_grad_desc_2 = m_hat_grad_desc_2.*(theta_star_hat_grad_desc_2(:,2) + a_m);

final_estimations_grad_desc_2 = [m_hat_grad_desc_2(end), k_hat_grad_desc_2(end), b_hat_grad_desc_2(end)];

% Create the graphical representations of the output x(t) (first component
% of the state vector x of the system) of the real system, the estimate x_hat(t)
% from the designed model, as well as the difference error_x(t) = x(t) - x_hat(t)
% which represents the modeling error of the output with the parallel topology
% model. Additionally, create the graphical representations of the real-time
% estimators m_hat(t), k_hat(t) and b_hat(t) for the time the simulation was
% run.

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t_grad_desc_2, x_grad_desc_2(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_grad_desc_2, x_hat_grad_desc_2(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_grad_desc_2, x_grad_desc_2(:,1)-x_hat_grad_desc_2(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Gradient Descent - Sinusoidal Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_grad_desc_2, m_hat_grad_desc_2, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_2, m_real*ones(size(t_grad_desc_2)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Gradient Descent - Sinusoidal Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_grad_desc_2, k_hat_grad_desc_2, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_2, k_real*ones(size(t_grad_desc_2)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Gradient Descent - Sinusoidal Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_grad_desc_2, b_hat_grad_desc_2, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_grad_desc_2, b_real*ones(size(t_grad_desc_2)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Gradient Descent - Sinusoidal Input', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on
