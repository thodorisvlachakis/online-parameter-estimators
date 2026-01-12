function dxdt = rollAngleDynamicSystem(t, x, a1, a2, a3, b, u, d)
% This function defines the state equations of the nonlinear dynamic system
% of the roll angle of an aircraft. The dynamic roll angle system is described
% by the following non-linear differential equation:
% r''(t) = - a1 * r'(t) - a2 * sin(r(t)) + a3 * (r'(t))^2 * sin(2*r(t)) + b * u(t) + d(t),
% where r(t) [rad] is the roll angle, a1,a2,a3 > 0 and b > 0 are constant
% parameters of the dynamic system, u(t) is the system control input and
% d(t) is external disturbances.
% The function rollAngleDynamicSystem defines the state equations that 
% describe the dynamic system of the roll angle of an aircraft, assuming the
% state variables x1(t) = r(t) and x2(t) = r'(t).
%
% The function takes as input the following:
% - t: the simulation time
% - x: the state vector
% - a1: a positive number which represents the parameter a1
% - a2: a positive number which represents the parameter a2
% - a3: a positive number which represents the parameter a3
% - b: a positive number which represents the parameter b
% - u: a function of time that represents the control input of the system
%      or a state feedback controller, which is actually a function of both
%      time t and state vector x.
% - d: a function of time or a vector representing the possible external
%      disturbances
    
    if nargin < 8
        d = 0;
    end

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

    if( isa(d, 'function_handle') )
        % if u is a function of time compute its value for the specific time t.
        d_val = d(t);
    else
        d_val = d;
    end

    dxdt = [x(2);
            -a1*x(2) - a2*sin(x(2)) + (a3 * x(2)^2 * sin(2*x(1))) + b*u_val + d_val];

end