% Check Exercise_2a_Roll_Angle_Sys_Controller_Implementation.m first.
% Tasks b and c below.

% ------------------------ Exercise 2 ------------------------
% In this exercise we consider the nonlinear dynamic system of the roll angle
% of an aircraft, which is described by the following nonlinear differential
% equation: 
% r''(t) = - a1 * r'(t) - a2 * sin(r(t)) + a3 * (r'(t))^2 * sin(2*r(t)) + b * u(t) + d(t),
% where r(t) [rad] is the roll angle, a1,a2,a3 > 0 and b > 0 are constant
% parameters of the dynamic system, u(t) is the system control input and
% d(t) is external disturbances.
% In this exercise, in task a, we will implement a state feedback controller
% in order to achieve a desired control objective for the roll angle r(t) of
% the aircraft system. In tasks b and c, we consider that the system parameters
% a1, a2, a3, b are unknown and we will design real-time estimators, using 
% composite topology model design and utilizing the Lyapunov method for the
% design of the estimators.
% In the following tasks, we will perform a simulation of the roll angle
% system using an ode solver in order to simulate the response of this system.
% The simulation of the system will make available the state vector, which
% we assume is available for measurement in the design of real-time parameter
% estimator. Obviously, simulating the system requires knowledge of it, i.e. 
% the actual values ​​of the parameters a1, a2, a3 and b. For our experiments,
% we will define the actual values ​​of the parameters and design real-time estimators 
% for them, using the Lyapunov method through composite topology (as we have
% nonlinear dynamic system).

addpath(fullfile(pwd, '..', 'functions'));

% Declare the actual values of the parameters a1, a2, a3 and b
a1 = 1.315;
a2 = 0.725;
a3 = 0.225;
b = 1.175;

% Define the reference trajectory for the roll angle
T_ramp = 5;
T_hold = 10;
rd_bar = pi / 10 ;

rd = desiredRollAngleTrajectory(rd_bar, T_ramp, T_hold);

% Declare the initial conditions for the roll angle system
% Rememeber that the state vector is x(t) = [r(t) r'(t)]^T.
x0 = [0; 0];

% Define the simulation time to be 20 sec and declare a strict condtion
% for the integration step Δt in order to achieve more accurate results.
tspan = [0 20];
integration_step_lim = 1e-3;

options = odeset('MaxStep',integration_step_lim);

figno = 0;

% b) In this task consider that the external disturbances are equal to zero. 
% We will design a real-time estimator of the unknown parameters a1, a2, a3
% and b of the real system, using the composite structure for model design
% and exploiting Lyapunov method to define real-time estimators.
d=0;
%d = @(t) 15 * sin(0.5*t) ;

% Define the controller that achieves the control objective for the roll 
% angle r(t), as it resulted from task a) in Exercise_2a_Roll_Angle_Sys_Controller_Implementation.m

k1 = 4;
k2 = 4;
lamda = 0.5;
rho = 3;
phi_inf = 0.01;
phi0 = 0.1;
controller = @(t, x) rollAngleSysController(t,x,rd,phi0,phi_inf,rho,k1,k2,lamda) ;

% Declare the initial conditions for the model and the real-time estimator
x0_hat = [0; 0];
theta0_hat = [2; 0.2; 0.6; 0.05];

% Declare the initial conditions for the extended state vector X which will
% be used to perform the simulation of the real system and the model and
% compute the real-time estimators.
X0 = [x0; x0_hat; theta0_hat];

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [a1 a2 a3 b]^T
Gamma = diag([44, 199, 3000, 38]);
theta_m = [0.3; 0.055];

[t, X] = ode45(@(t,X) composite_structure_model_adaptive_estimation_RollAngleSys(t,X,a1,a2,a3,b,controller,Gamma,theta_m,d), tspan, X0, options);

x = X(:, 1:2);
x_hat = X(:, 3:4);
theta_hat = X(:, 5:8);

final = theta_hat(end, :);

% Create the graphical representations of the roll angle r(t) (first component
% of the state vector x of the system) of the real system, the estimate r_hat(t)
% from the designed model, as well as the difference error_r(t) = r(t) - r_hat(t)
% Additionally, create the graphical representations of the real-time
% estimators a1_hat(t), a2_hat(t), a3_hat(t) and b_hat(t) for the time the
% simulation was run.

figno = figno + 1;

figure(figno)
clf
plot(t, x(:,1), 'DisplayName','r(t)', 'Color', 'black')
hold on
plot(t, x_hat(:,1), 'DisplayName','$\hat{r}$ (t)', 'Color', 'green')
hold on
plot(t, x(:,1)-x_hat(:,1), 'DisplayName', '$e_r(t) = r(t) - \hat{r}(t)$', 'Color', 'red')
title('r(t), $\hat{r}(t)$ and $e_r(t) = r(t) - \hat{r}(t)$ , Composite Topology-Lyapunov Method, d(t)=0', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('r(t), $\hat{r}(t), e_r(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t, theta_hat(:,1), 'DisplayName','$\hat{a}_1$ (t)', 'Color', 'blue')
hold on
plot(t, a1*ones(size(t)), 'DisplayName','$a_1$', 'Color', 'red')
hold on
title('Graph of $\hat{a}_1(t)$, Composite Topology-Lyapunov Method, d(t)=0', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_1$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t, theta_hat(:,2), 'DisplayName','$\hat{a}_2$ (t)', 'Color', 'blue')
hold on
plot(t, a2*ones(size(t)), 'DisplayName','$a_2$', 'Color', 'red')
title('Graph of $\hat{a}_2(t)$, Composite Topology-Lyapunov Method, d(t)=0', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_2$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t, theta_hat(:,3), 'DisplayName','$\hat{a}_3$ (t)', 'Color', 'blue')
hold on
plot(t, a3*ones(size(t)), 'DisplayName','$a_3$', 'Color', 'red')
title('Graph of $\hat{a}_3(t)$, Composite Topology-Lyapunov Method, d(t)=0', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_3$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t, theta_hat(:,4), 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t, b*ones(size(t)), 'DisplayName','b', 'Color', 'red')
title('Graph of $\hat{b}(t)$, Composite Topology-Lyapunov Method, d(t)=0', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

% b) In this task, consider the external disturbances as: d(t) = 15 * sin(0.5*t)
% and re-implement the real-time estimator of task b) to study the effect
% of the disturbances.
d = @(t) 15 * sin(0.5*t) ;

% Define the controller that achieves the control objective for the roll 
% angle r(t), as it resulted from task a) in Exercise_2a_Roll_Angle_Sys_Controller_Implementation.m

controller = @(t, x) rollAngleSysController(t,x,rd,phi0,phi_inf,rho,k1,k2,lamda) ;

% Declare the initial conditions for the model and the real-time estimator
x0_hat = [0; 0];
theta0_hat = [2; 0.2; 0.6; 0.05];

% Declare the initial conditions for the extended state vector X which will
% be used to perform the simulation of the real system and the model and
% compute the real-time estimators.
X0 = [x0; x0_hat; theta0_hat];

% Define gain (learning rate) for the estimator of each parameter of the
% vector theta = [a1 a2 a3 b]^T
Gamma = diag([44, 199, 3000, 38]);
theta_m = [0.3; 0.055];

[t_disturbed, X_disturbed] = ode45(@(t,X) composite_structure_model_adaptive_estimation_RollAngleSys(t,X,a1,a2,a3,b,controller,Gamma,theta_m,d), tspan, X0, options);

x_disturbed = X_disturbed(:, 1:2);
x_hat_disturbed = X_disturbed(:, 3:4);
theta_hat_disturbed = X_disturbed(:, 5:8);

final_disturbed = theta_hat_disturbed(end, :);

% Create the graphical representations of the real-time estimators a1_hat(t),
% a2_hat(t), a3_hat(t) and b_hat(t) as well as graphical representations of
% the absolute errors |a1_hat(t) - a1|, |a2_hat(t) - a2|, |a3_hat(t) - a3|
% and |b_hat(t) - b| , for the time the simulation was run.

figno = figno + 1;

figure(figno)
clf
subplot(2,2,1)
plot(t_disturbed, theta_hat_disturbed(:,1), 'DisplayName','$\hat{a}_1$ (t)', 'Color', 'blue')
hold on
plot(t_disturbed, a1*ones(size(t_disturbed)), 'DisplayName','$a_1$', 'Color', 'red')
hold on
title('Graph of $\hat{a}_1(t)$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_1$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_disturbed, abs(theta_hat_disturbed(:,1) - a1*ones(size(t_disturbed))), 'DisplayName', '$|\hat{a}_1(t) - a_1|$', 'Color',[0.6350 0.0780 0.1840])
hold on
title('$|\hat{a}_1(t) - a_1|$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$|\hat{a}_1(t) - a_1|$', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_disturbed, theta_hat_disturbed(:,2), 'DisplayName','$\hat{a}_2$ (t)', 'Color', 'blue')
hold on
plot(t_disturbed, a2*ones(size(t_disturbed)), 'DisplayName','$a_2$', 'Color', 'red')
title('Graph of $\hat{a}_2(t)$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_2$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_disturbed, abs(theta_hat_disturbed(:,2) - a2*ones(size(t_disturbed))), 'DisplayName', '$|\hat{a}_2(t) - a_2|$', 'Color',[0.6350 0.0780 0.1840])
hold on
title('$|\hat{a}_2(t) - a_2|$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$|\hat{a}_2(t) - a_2|$', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

figno = figno + 1;
figure(figno)
clf
subplot(2,2,1)
plot(t_disturbed, theta_hat_disturbed(:,3), 'DisplayName','$\hat{a}_3$ (t)', 'Color', 'blue')
hold on
plot(t_disturbed, a3*ones(size(t_disturbed)), 'DisplayName','$a_3$', 'Color', 'red')
title('Graph of $\hat{a}_3(t)$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{a}_3$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,2)
plot(t_disturbed, abs(theta_hat_disturbed(:,3) - a3*ones(size(t_disturbed))), 'DisplayName', '$|\hat{a}_3(t) - a_3|$', 'Color',[0.6350 0.0780 0.1840])
hold on
title('$|\hat{a}_3(t) - a_3|$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$|\hat{a}_3(t) - a_3|$', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,3)
plot(t_disturbed, theta_hat_disturbed(:,4), 'DisplayName','$\hat{b}$ (t)', 'Color', 'blue')
hold on
plot(t_disturbed, b*ones(size(t_disturbed)), 'DisplayName','b', 'Color', 'red')
title('Graph of $\hat{b}(t)$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$\hat{b}$ (t)', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on

subplot(2,2,4)
plot(t_disturbed, abs(theta_hat_disturbed(:,4) - b*ones(size(t_disturbed))), 'DisplayName', '$|\hat{b}(t) - b|$', 'Color',[0.6350 0.0780 0.1840])
hold on
title('$|\hat{b}(t) - b|$, Composite Topology-Lyapunov Method, d(t)=0.15$\cdot$sin(0.5$\cdot$t)', 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
ylabel('$|\hat{b}(t) - b|$', 'Interpreter','Latex', 'FontWeight','bold')
legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
grid on
