function rd = desiredRollAngleTrajectory(rd_bar, T_ramp, T_hold)
% This function creates a smooth reference trajectory rd(t) which specifies
% a control objective (goal) for the roll angle r(t) of an aircraft, at a
% certain time interval. The control objective is to drive the roll angle
% r(t) from the initial value r(0) = 0 to a desired (target) value rd_bar 
% and then return to zero angle. This trajectory (rd(t): 0 -> rd_bar -> 0)
% that we want the roll angle to follow is implemented through a smooth 
% reference trajectory rd(t) that satisfies this objective, in a certain
% amount of time (e.g. 20 sec).
% The function desiredRollAngleTrajectory this smooth trajectory rd(t) 
% considering a time interval of "rising" the aircraft's roll angle from
% zero to the target value, a time interval of "staying" the roll angle at
% this value and finally a time interval of "descending" the roll angle back
% to zero. The function, then, constructs this smooth trajectory by creating
% a smooth branch function of time rd(t), which is modeled by a 5th degree
% polynomial in the "rising" time interval, a constant function with value
% rd_bar in the "staying" time interval, and a 5th degree polynomial in the
% "descending" time interval. The use of 5th degree polynomial functions is
% beneficial in creating a smooth transition (rd(t): 0 -> rd_bar -> 0) in a
% controlled manner.
% The function takes as input the following:
% - rd_bar: a positive number which represents the target value for the
%           roll angle.
% - T_ramp: a positive number which represents both the "rising" and the
%           "descending" time interval.
% - T_hold: a positive number which represnts the "staying" time interval.
%
% Outputs: Returns a function handle that represents the desired trajectory

    % Define the 5th degree polynomial function to be used for the transitions
    % (rd(t): 0 -> rd_bar and rd(t): rd_bar -> 0)

    % quintic_polynomial = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5
    % Compute the parameters of the quintic polynomial, by solving a system
    % of differential equations, which results from the initial conditions
    % that the polynomial function must satisfy as well as its first and
    % second derivatives.

    % Initial conditions: rd(0)=0, rd(T_ramp)=rd_bar, rd'(0)=0, rd'(T_ramp)=0,
    % rd''(0)=0, rd''(T_ramp)=0. Construct the system of differential equations.

    A = [1 0 0 0 0 0;
         1, T_ramp, T_ramp^2, T_ramp^3, T_ramp^4, T_ramp^5;
         0 1 0 0 0 0;
         0, 1, 2*T_ramp, 3*T_ramp^2, 4*T_ramp^3, 5*T_ramp^4;
         0 0 2 0 0 0;
         0, 0, 2, 6*T_ramp, 12*T_ramp^2, 20*T_ramp^3];

    b = [0; rd_bar; 0; 0; 0; 0];
    
    quintic_pol_coeff = A\b ;

    quintic_polynomial = @(t) quintic_pol_coeff(1) + quintic_pol_coeff(2)*t + quintic_pol_coeff(3)*t.^2 + ...
                   quintic_pol_coeff(4)*t.^3 + quintic_pol_coeff(5)*t.^4 + quintic_pol_coeff(6)*t.^5;

    rd = @(t) (t >= 0 & t <= T_ramp) .* quintic_polynomial(t) + ... 
              (t > T_ramp & t <= T_ramp + T_hold) .* rd_bar + ... 
              (t > T_ramp + T_hold & t <= 2*T_ramp + T_hold) .* quintic_polynomial(2*T_ramp + T_hold - t) + ...
              (t > 2*T_ramp + T_hold) .* 0;

end