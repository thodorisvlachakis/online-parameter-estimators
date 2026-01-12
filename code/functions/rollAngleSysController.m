function [u, phi, a] = rollAngleSysController(t, x, rd, phi0, phi_inf, rho, k1, k2, lamda)
% This function implements a state feedback controller u(t) = u(r(t), r'(t))
% to achieve the controll objective (goal) for the roll angle system, assuming
% that the external disturbances d(t) is equal to 0.
% The dynamic roll angle system is described by the following non-linear
% differential equation:
% r''(t) = - a1 * r'(t) - a2 * sin(r(t)) + a3 * (r'(t))^2 * sin(2*r(t)) + b * u(t) + d(t),
% where r(t) [rad] is the roll angle, a1,a2,a3 > 0 and b > 0 are constant
% parameters of the dynamic system, u(t) is the system control input and
% d(t) is external disturbances.
% The control objective is to drive the roll angle r(t) from the initial value
% r(0) = 0 to a desired (target) value rd_bar and then return to zero angle.
% This trajectory (rd(t): 0 -> rd_bar -> 0) that we want the roll angle to
% follow is implemented through a smooth reference trajectory rd(t) that
% satisfies this objective, in a certain amount of time (e.g. 20 sec). So,
% the state feedback controller u(t) = u(r(t), r'(t)) drives the system to
% follow the reference trajectory rd(t), which is defined in advance so as
% to satisfy the desired specifications. The state feedback controller
% u(t) = u(r(t), r'(t)) the function implements is given by the following
% formula:
% z1(t) = (r(t) - rd(t)) / phi(t) , a(t) = -k1 * T(z1(t)),
% z2(t) = (r'(t) - a(t)) / rho , u(t) = -k2 * T(z2(t)),
% where k1,k2 > 0 (gains), phi(t) = (phi0 - phi_inf)*exp(-lamda*t) + phi_inf
% with parameters phi0 > phi_inf > 0 , lamda > 0 , phi0 >> |r(0) - rd(0)|,
% T(z) = ln((1+z)/(1-z)) and rho >> |r'(0) - a0|
%
% The function takes as input the following:
% - t: the simulation time.
% - x: the state vector, whose first component represent the roll angle
%      r(t) and the second represents the derivative r'(t).
% - rd: a function of time that represents the desired trajectory for r(t).
% - phi0: a positive number which a parameter of phi(t) function.
% - phi_inf: a positive number which a parameter of phi(t) function.
% - rho: a positive number which represents a parameter of the controller.
% - k1: a positive number which represents a gain.
% - k2: a positive number which represents a gain.

    roll_angle = x(1);
    roll_angle_dot = x(2);

    if( isa(rd, 'function_handle') )
        % if rd is a function of time compute its value for the specific time t.
        rd_val = rd(t);
    else
        rd_val = rd;
    end
    
    %{
    T = @(z) log( ( 1+z ) / ( 1-z ) );

    phi = @(t) (phi0 - phi_inf) * exp(-lamda*t) + phi_inf;

    z1 = @(t) (r(t) - rd(t)) / phi(t);
    a = @(t) -k1 * T(z1(t)) ;
    %}

    T = @(z) log( ( 1+z ) / ( 1-z ) );

    phi = (phi0 - phi_inf) * exp(-lamda*t) + phi_inf;

    z1 = (roll_angle - rd_val) / phi;
    a = -k1 * T(z1) ;
    
    % Check that the selected values for the parameters of the controller
    % satisfy the necessary conditions.
    if t==0
        if(rho < 10 * abs(roll_angle - a))
            print("The value of the 'rho' parameter of the controller is not as large as needed.")
            return
        end
    end

    z2 = (roll_angle_dot - a) / rho;
    u = -k2 * T(z2);

end