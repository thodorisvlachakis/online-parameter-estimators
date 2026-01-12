% Task a) below.

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

% a) In this task, we will implement a state feedback controller u(t) = u(r(t), r'(t))
% to achieve the controll objective (goal) for the roll angle system when
% external disturbances are considered zero (d(t)=0). The control objective
% is to drive the roll angle r(t) from the initial value r(0) = 0 to a desired
% (target) value rd_bar and then return to zero angle. This trajectory
% (rd(t): 0 -> rd_bar -> 0) that we want the roll angle to follow is implemented
% through a smooth reference trajectory rd(t) that satisfies this objective,
% in a certain amount of time (e.g. 20 sec).
% We will use an ode solver in order to simulate the system for a simulation
% time period of 20 sec, with an appropriate integration step for accurate
% results, considering the appropriate state feedback controller, which we
% implement so that it achieves the control objective, as input u = u(t) of
% the roll angle system and we will obtain the response of the closed-loop
% system.
 
% Define the parameters of the controller to be implemented by the function
% rollAngleSysController.
% Try out different values for phi_inf.
k1 = 4;
k2 = 4;
lamda = 0.5;
rho = 3;
phi_inf = [1; 0.4; 0.1; 0.01] + abs(x0(1)-rd(0));
phi0 = 10*phi_inf;

for i=1:length(phi_inf)
    controller = @(t,x) rollAngleSysController(t,x,rd,phi0(i),phi_inf(i),rho,k1,k2,lamda);

    % Perform simulation of the closed-loop system, applying the controller
    % that achieves the controll objective.
 
    [t, x] = ode45(@(t,x) rollAngleDynamicSystem(t, x, a1, a2, a3, b, controller), tspan, x0, options);

    u = zeros(size(t));
    phi = zeros(size(t));
    a = zeros(size(t));

    for j=1:length(t)
        [u(j), phi(j), a(j)] = rollAngleSysController(t(j), x(j,:), rd, phi0(i), phi_inf(i), rho, k1, k2, lamda);
    end

    % Create the graphical representations of the roll angle r(t) (first component
    % of the state vector x of the system) of the real system along with the
    % desired trajectory rd(t).
    % Also create graphical representations of the signals |rd(t) - r(t)|
    % and phi(t) in one graph as well as graphical representation of the signal
    % |r'(t) - a(t)| along with the constant value of rho in an another graph
    % in order to reach some conclusions.
    figno = figno + 1;
    
    figure(figno)
    clf
    subplot(3,1,1)
    plot(t, x(:,1), 'DisplayName','r(t)', 'Color', 'blue')
    hold on
    plot(t, rd(t), 'DisplayName','$r_d$ (t)', 'Color', 'red')
    hold on
    title(['r(t) and $r_d(t)$, Roll Angle Controll ($\phi_\infty$ = {',num2str(phi_inf(i)),'})'], 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
    xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
    ylabel('r(t), $r_d(t)$ [rad]', 'Interpreter','Latex', 'FontWeight','bold')
    legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
    grid on
    hold off
    
    subplot(3,1,2)
    plot(t, phi, 'DisplayName','$\phi(t)$', 'Color', 'blue')
    hold on
    plot(t, abs(x(:,1) - rd(t)), 'DisplayName','$|r(t) - r_d(t)|$', 'Color', 'red')
    title(['$\phi(t)$ and $|r(t) - r_d(t)|$, State Feedback Controller ($\phi_\infty$ = {',num2str(phi_inf(i)),'})'], 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
    xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
    ylabel('$\phi(t)$, $|r(t) - r_d(t)|$', 'Interpreter','Latex', 'FontWeight','bold')
    legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
    grid on
    
    subplot(3,1,3)
    plot(t, rho*ones(size(t)), 'DisplayName','$\rho$', 'Color', 'blue')
    hold on
    plot(t, abs(x(:,2) - a), 'DisplayName','$|\dot{r}(t) - a(t)|$', 'Color', 'red')
    title(['$|\dot{r}(t) - a(t)|$, State Feedback Controller ($\phi_\infty$ = {',num2str(phi_inf(i)),'})'], 'Interpreter','Latex', 'FontWeight', 'bold', 'FontSize', 12)
    xlabel('t [sec]','Interpreter','Latex', 'FontWeight', 'bold')
    ylabel('$|\dot{r}(t) - a(t)|$', 'Interpreter','Latex', 'FontWeight','bold')
    legend('Interpreter','Latex','FontWeight','bold', 'Location', 'Best')
    grid on

end

% Continue with the smallest phi0 and phi_inf values, since the accuracy of
% tracking the desired trajectory is improved.
phi0 = phi0(4);
phi_inf = phi_inf(4);