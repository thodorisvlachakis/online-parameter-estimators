% Tasks b and c below

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

% b) In this task we will design a real-time estimator of the unknown parameters
% m, k and b of the real system, using initially the parallel structure (topology)
% for model design and exploiting Lyapunov method to define real-time
% estimators and then the composite structure for model design and exploiting
% Lyapunov method to define real-time estimators.

% Define the input of the system
u = @(t) 2.5*sin(t) ;

% Declare the initial conditions for the model and the real-time estimator
x0_hat = [0; 0];
theta0_hat = [2; 0.2; 0.5];

% Declare the initial conditions for the extended state vector X which will
% be used to perform the simulation of the real system and the model and
% compute the real-time estimators.
X0 = [x0; x0_hat; theta0_hat];

% b1) Design of real-time estimator of the unknown parameters m, k, b of the
% system, using a parallel structure model and the Lyapunov method for the
% design.

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [k/m b/m 1/m]^T
Gamma = diag([15, 1.5, 7]);

[t_parallel_topology, X_parallel_topology] = ode45(@(t,X) parallel_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma), tspan, X0, options);

x_parallel_topology = X_parallel_topology(:, 1:2);
x_hat_parallel_topology = X_parallel_topology(:, 3:4);
theta_hat_parallel_topology = X_parallel_topology(:, 5:7);

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
% and solving the equations resulting from the corresponding components.
m_hat_parallel_topology = 1./theta_hat_parallel_topology(:,3);
k_hat_parallel_topology = m_hat_parallel_topology.*theta_hat_parallel_topology(:,1);
b_hat_parallel_topology = m_hat_parallel_topology.*theta_hat_parallel_topology(:,2);

final_estimations_par_top = [m_hat_parallel_topology(end), k_hat_parallel_topology(end), b_hat_parallel_topology(end)];

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
plot(t_parallel_topology, x_parallel_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_parallel_topology, x_hat_parallel_topology(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_parallel_topology, x_parallel_topology(:,1)-x_hat_parallel_topology(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('Graphical Representations of x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_parallel_topology, m_hat_parallel_topology, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, m_real*ones(size(t_parallel_topology)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_parallel_topology, k_hat_parallel_topology, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, k_real*ones(size(t_parallel_topology)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_parallel_topology, b_hat_parallel_topology, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, b_real*ones(size(t_parallel_topology)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on


% b2) Design of real-time estimator of the unknown parameters m, k, b of the
% system, using a composite model structure and the Lyapunov method for the
% design.

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [k/m b/m 1/m]^T
Gamma = diag([15, 1.5, 7]);
theta_m = [2; 2];

[t_composite_topology, X_composite_topology] = ode45(@(t,X) composite_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma,theta_m), tspan, X0, options);

x_composite_topology = X_composite_topology(:, 1:2);
x_hat_composite_topology = X_composite_topology(:, 3:4);
theta_hat_composite_topology = X_composite_topology(:, 5:7);

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
% and solving the equations resulting from the corresponding components.
m_hat_composite_topology = 1./theta_hat_composite_topology(:,3);
k_hat_composite_topology = m_hat_composite_topology.*theta_hat_composite_topology(:,1);
b_hat_composite_topology = m_hat_composite_topology.*theta_hat_composite_topology(:,2);

final_estimations = [m_hat_composite_topology(end), k_hat_composite_topology(end), b_hat_composite_topology(end)];

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
plot(t_composite_topology, x_composite_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_composite_topology, x_composite_topology(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_composite_topology, x_composite_topology(:,1)-x_hat_composite_topology(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('Graphical Representations of x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_composite_topology, m_hat_composite_topology, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, m_real*ones(size(t_composite_topology)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_composite_topology, k_hat_composite_topology, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, k_real*ones(size(t_composite_topology)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_composite_topology, b_hat_composite_topology, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, b_real*ones(size(t_composite_topology)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on


% c) In this task, we will repeat the process of designing the real-time estimators
% of the unknown parameters m, k and b of the system in case of both a parallel
% structure and a composite structure model design, exactly as in task b), 
% but now we consider that the output x(t) (which is the first component
% of the state vector x of the system) of the real system is measured with
% noise n(t) = n0 * sin(2*pi*f0*t), where n0=0.25 and f0=20.

% Define the input of the system
u = @(t) 2.5*sin(t) ;

% Define the noise that is measured in measurements of the output x(t).
n0 = 0.25;
f0 = 20;
n = @(t) n0*sin(2*pi*f0*t) ;

% c1) Design of real-time estimator of the unknown parameters m, k, b of the
% system, using a parallel structure model and the Lyapunov method for the
% design.

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [k/m b/m 1/m]^T
Gamma = diag([15, 1.5, 7]);

[t_parallel_topology_noisy, X_parallel_topology_noisy] = ode45(@(t,X) parallel_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma,n), tspan, X0, options);

x_parallel_topology_noisy = X_parallel_topology_noisy(:, 1:2);
x_hat_parallel_topology_noisy = X_parallel_topology_noisy(:, 3:4);
theta_hat_parallel_topology_noisy = X_parallel_topology_noisy(:, 5:7);

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
% and solving the equations resulting from the corresponding components.
m_hat_parallel_topology_noisy = 1./theta_hat_parallel_topology_noisy(:,3);
k_hat_parallel_topology_noisy = m_hat_parallel_topology_noisy.*theta_hat_parallel_topology_noisy(:,1);
b_hat_parallel_topology_noisy = m_hat_parallel_topology_noisy.*theta_hat_parallel_topology_noisy(:,2);

final_estimations_noisy_par_top = [m_hat_parallel_topology_noisy(end), k_hat_parallel_topology_noisy(end), b_hat_parallel_topology_noisy(end)];

% Create the graphical representations of the output x(t) (first component
% of the state vector x of the system) of the real system, the estimate x_hat(t)
% from the designed model, as well as the difference error_x(t) = x(t) - x_hat(t)
% which represents the modeling error of the output with the parallel topology
% model. Additionally, create the graphical representations of the real-time
% estimators m_hat(t), k_hat(t) and b_hat(t) for the time the simulation was
% run. Compare the results with the corresponding ones without noise.

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t_parallel_topology, x_parallel_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_parallel_topology, x_hat_parallel_topology(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_parallel_topology, x_parallel_topology(:,1)-x_hat_parallel_topology(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_parallel_topology, x_parallel_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_parallel_topology_noisy, x_hat_parallel_topology_noisy(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_parallel_topology_noisy, x_parallel_topology(:,1)-x_hat_parallel_topology_noisy(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Parallel Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_parallel_topology, m_hat_parallel_topology, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, m_real*ones(size(t_parallel_topology)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_parallel_topology_noisy, m_hat_parallel_topology_noisy, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology_noisy, m_real*ones(size(t_parallel_topology_noisy)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Parallel Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

figno = figno + 1;
figure(figno)
clf
subplot(2,2,1)
plot(t_parallel_topology, k_hat_parallel_topology, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, k_real*ones(size(t_parallel_topology)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_parallel_topology_noisy, k_hat_parallel_topology_noisy, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology_noisy, k_real*ones(size(t_parallel_topology_noisy)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Parallel Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_parallel_topology, b_hat_parallel_topology, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology, b_real*ones(size(t_parallel_topology)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Parallel Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_parallel_topology_noisy, b_hat_parallel_topology_noisy, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_parallel_topology_noisy, b_real*ones(size(t_parallel_topology_noisy)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Parallel Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on


% c2) Design of real-time estimator of the unknown parameters m, k, b of the
% system, using a composite model structure and the Lyapunov method for the
% design.

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [k/m b/m 1/m]^T
Gamma = diag([15, 1.5, 7]);
theta_m = [2; 2];

[t_composite_topology_noisy, X_composite_topology_noisy] = ode45(@(t,X) composite_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma,theta_m,n), tspan, X0, options);

x_composite_topology_noisy = X_composite_topology_noisy(:, 1:2);
x_hat_composite_topology_noisy = X_composite_topology_noisy(:, 3:4);
theta_hat_composite_topology_noisy = X_composite_topology_noisy(:, 5:7);

% Now, compute the real-time estimators for the parameters m, k and b
% considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
% and solving the equations resulting from the corresponding components.
m_hat_composite_topology_noisy = 1./theta_hat_composite_topology_noisy(:,3);
k_hat_composite_topology_noisy = m_hat_composite_topology_noisy.*theta_hat_composite_topology_noisy(:,1);
b_hat_composite_topology_noisy = m_hat_composite_topology_noisy.*theta_hat_composite_topology_noisy(:,2);

final_estimations_noisy = [m_hat_composite_topology_noisy(end), k_hat_composite_topology_noisy(end), b_hat_composite_topology_noisy(end)];

% Create the graphical representations of the output x(t) (first component
% of the state vector x of the system) of the real system, the estimate x_hat(t)
% from the designed model, as well as the difference error_x(t) = x(t) - x_hat(t)
% which represents the modeling error of the output with the parallel topology
% model. Additionally, create the graphical representations of the real-time
% estimators m_hat(t), k_hat(t) and b_hat(t) for the time the simulation was
% run. Compare the results with the corresponding ones without noise.

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t_composite_topology, x_composite_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_composite_topology, x_hat_composite_topology(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_composite_topology, x_composite_topology(:,1)-x_hat_composite_topology(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_composite_topology, x_composite_topology(:,1), 'DisplayName','x(t)', 'Color', 'black')
hold on
plot(t_composite_topology_noisy, x_hat_composite_topology_noisy(:,1), 'DisplayName','$\hat{x}$ (t)', 'Color', 'green')
hold on
plot(t_composite_topology_noisy, x_composite_topology(:,1)-x_hat_composite_topology_noisy(:,1), 'DisplayName', '$e_x(t) = x(t) - \hat{x}(t)$', 'Color', 'red')
title('x(t), $\hat{x}(t)$ and $e_x(t) = x(t) - \hat{x}(t)$ , Composite Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('x(t), $\hat{x}(t), e_x(t)$ [length units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_composite_topology, m_hat_composite_topology, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, m_real*ones(size(t_composite_topology)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_composite_topology_noisy, m_hat_composite_topology_noisy, 'DisplayName','$\hat{m}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology_noisy, m_real*ones(size(t_composite_topology_noisy)), 'DisplayName','m', 'Color', 'red')
title('Graphical Representation of $\hat{m}(t)$, Composite Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{m}$ (t) [mass units]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

figno = figno + 1;
figure(figno)
clf
subplot(2,2,1)
plot(t_composite_topology, k_hat_composite_topology, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, k_real*ones(size(t_composite_topology)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_composite_topology_noisy, k_hat_composite_topology_noisy, 'DisplayName','$\hat{k}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology_noisy, k_real*ones(size(t_composite_topology_noisy)), 'DisplayName','k', 'Color', 'red')
title('Graphical Representation of $\hat{k}$(t), Composite Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{k}$ (t) [N/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_composite_topology, b_hat_composite_topology, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology, b_real*ones(size(t_composite_topology)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Composite Topology', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_composite_topology_noisy, b_hat_composite_topology_noisy, 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_composite_topology_noisy, b_real*ones(size(t_composite_topology_noisy)), 'DisplayName','b', 'Color', 'red')
title('Graphical Representation of $\hat{b}(t)$, Composite Topology - Sinusoidal Noise Effect', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t) [N $\cdot$ s/m]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on


% Final Step: He studied the effect of varying the noise amplitude n0 on the
% accuracy of the estimated parameters with real-time estimators of both the
% parallel structure and the composite structure with the Lyapunov method, 
% at the end of the simulation time, when theoretically the estimators are
% closest to converging to the true value of the corresponding parameter.

n0_values = linspace(0.1,50,20);
absolute_estimation_errors_par_top = zeros(length(n0_values), 3);
euclidean_estimation_error_par_top = zeros(length(n0_values), 1);
absolute_estimation_errors_com_top = zeros(length(n0_values), 3);
euclidean_estimation_error_com_top = zeros(length(n0_values), 1);

for i=1:length(n0_values)
    n = @(t) i*sin(2*pi*f0*t) ;

    % c1) Design of real-time estimator of the unknown parameters m, k, b of the
    % system, using a parallel structure model and the Lyapunov method for the
    % design.
    
    [t_par_top_noisy, X_par_top_noisy] = ode45(@(t,X) parallel_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma,n), tspan, X0, options);
    
    x_par_top_noisy = X_par_top_noisy(:, 1:2);
    x_hat_par_top_noisy = X_par_top_noisy(:, 3:4);
    theta_hat_par_top_noisy = X_par_top_noisy(:, 5:7);
    
    % Now, compute the real-time estimators for the parameters m, k and b
    % considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
    % and solving the equations resulting from the corresponding components.
    m_hat_par_top_noisy = 1./theta_hat_par_top_noisy(:,3);
    k_hat_par_top_noisy = m_hat_par_top_noisy.*theta_hat_par_top_noisy(:,1);
    b_hat_par_top_noisy = m_hat_par_top_noisy.*theta_hat_par_top_noisy(:,2);
    
    sys_params_par_top_noisy = [m_hat_par_top_noisy(end), k_hat_par_top_noisy(end), b_hat_par_top_noisy(end)];
    
    absolute_estimation_errors_par_top(i,:) = abs([m_real, k_real, b_real] - sys_params_par_top_noisy);
    euclidean_estimation_error_par_top(i) = norm([m_real, k_real, b_real] - sys_params_par_top_noisy);

    % c2) Design of real-time estimator of the unknown parameters m, k, b of the
    % system, using a composite model structure and the Lyapunov method for the
    % design.
    
    [t_com_top_noisy, X_com_top_noisy] = ode45(@(t,X) composite_structure_model_adaptive_estimation_MSD(t,X,m_real, b_real,k_real,u,Gamma,theta_m,n), tspan, X0, options);
    
    x_com_top_noisy = X_com_top_noisy(:, 1:2);
    x_hat_com_top_noisy = X_com_top_noisy(:, 3:4);
    theta_hat_com_top_noisy = X_com_top_noisy(:, 5:7);
    
    % Now, compute the real-time estimators for the parameters m, k and b
    % considering that the vector theta_hat is the estimation of theta = [k/m b/m 1/m]^T
    % and solving the equations resulting from the corresponding components.
    m_hat_com_top_noisy = 1./theta_hat_com_top_noisy(:,3);
    k_hat_com_top_noisy = m_hat_com_top_noisy.*theta_hat_com_top_noisy(:,1);
    b_hat_com_top_noisy = m_hat_com_top_noisy.*theta_hat_com_top_noisy(:,2);

    sys_params_com_top_noisy = [m_hat_com_top_noisy(end), k_hat_com_top_noisy(end), b_hat_com_top_noisy(end)];
    
    absolute_estimation_errors_com_top(i,:) = abs([m_real, k_real, b_real] - sys_params_com_top_noisy);
    euclidean_estimation_error_com_top(i) = norm([m_real, k_real, b_real] - sys_params_com_top_noisy);

end


% Graphical respresentations for parallel topology

figno = figno + 1;
figure(figno)
clf
subplot(2,2,1)
scatter(n0_values, absolute_estimation_errors_par_top(:,1), 'r', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $m$, Parallel Topology', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

subplot(2,2,2)
scatter(n0_values, absolute_estimation_errors_par_top(:,2), 'g', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $k$, Parallel Topology', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,3)
scatter(n0_values, absolute_estimation_errors_par_top(:,3), 'b', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex', 'FontWeight','bold')
ylabel('$|b - \hat{b}|$', 'Interpreter', 'latex', 'FontWeight','bold')
title('Absolute Error in Estimation of $b$, Parallel Topology', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,4)
scatter(n0_values, euclidean_estimation_error_par_top, 'k', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('Euclidean Error', 'Interpreter', 'latex','FontWeight','bold')
title('Total (Euclidean) Parameter Estimation Error, Parallel Topology', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

% Graphical representations for composite topology

figno = figno + 1;
figure(figno)
clf
subplot(2,2,1)
scatter(n0_values, absolute_estimation_errors_com_top(:,1), 'r', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $m$, Composite Topology', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)

subplot(2,2,2)
scatter(n0_values, absolute_estimation_errors_com_top(:,2), 'g', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex','FontWeight','bold')
title('Absolute Error in Estimation of $k$, Composite Topology', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,3)
scatter(n0_values, absolute_estimation_errors_com_top(:,3), 'b', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex', 'FontWeight','bold')
ylabel('$|b - \hat{b}|$', 'Interpreter', 'latex', 'FontWeight','bold')
title('Absolute Error in Estimation of $b$, Composite Topology', 'Interpreter', 'latex', 'FontWeight','bold','FontSize',12)

subplot(2,2,4)
scatter(n0_values, euclidean_estimation_error_com_top, 'k', 'filled')
hold on
grid on
xlabel('Noise Amplitude $\eta_0$', 'Interpreter', 'latex','FontWeight','bold')
ylabel('Euclidean Error', 'Interpreter', 'latex','FontWeight','bold')
title('Total (Euclidean) Parameter Estimation Error, Composite Topology', 'Interpreter', 'latex','FontWeight','bold','FontSize',12)
