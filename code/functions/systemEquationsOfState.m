function dxdt = systemEquationsOfState(t, x, A, B, u)
% This function defines the state equations of a dynamic system in the
% linear state-space representation: x' = A*x + B*u ,
% where x is the state vector and u is the control input.
% The function takes as input the following:
% - A 2D (n×n) array A which is the matrix A in linear state-space representation
% - A 2D (n×m) array B which is the matrix B in linear state-space representation
% - A vector x of length n which represents the state vector x
% - A vector u of length m which represents the control input. It can be
% scalar and in this case the array B is a vector of length n.
    
    if( isa(u, 'function_handle') )
        % if u is a function of time compute its value for the specific time t.
        u_val = u(t);
    else
        u_val = u;
    end
    
    dxdt = A * x + B * u_val ;

end