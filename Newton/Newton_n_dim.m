function [X] = newton_n_dim(tolerance_rss,initial_estimate,sym_variables,sym_equations)
%% newton method for solving a system of >=n nonlinear equations for n variables

% Given n equations, the function performs the newton method, converging
% to the exact solution.
% Given >n equations, the function converges to the solution which minimizes
% the least squared error of the given equations.
%
%input:     tolerance_norm :        maximum tolerable RSS of errors in output vector
%           initial estimate :      row vector of initial estimate
%           sym_variables:          row vector of n symbolic variables
%           sym_equations:          column vector of >=n symbolic equations
%
%output:    solution:               row vector of solution.
%
%assumptions:   1.  Input sym functions are differentiable
%               2.  Convergence is dependent on the functions.
%                  -check convergence constraints.
%                    -http://en.wikipedia.org/wiki/Newton's_method

%%   Example:
%        syms a b
%        F1 = a-15;             %(15 = a)
%        F2 = b^2-10;           %(10 = b^2)
%        tolerance = .01;
%        initial_est = [10,1];
%with n equations and n unknowns:
%        solution = newton_n_dim(tolerance,initial_est,[a,b],[F1;F2]);
%with >n equations and n unknowns:
%        F3 = sqrt(a^2 + b^2)-15.5; (third equation, (15.5 = sqrt(a^2 + b^2)))
%        solution = newton_n_dim(tolerance,initial_est,[a,b],[F1;F2;F3]);

%Kyle J. Drerup
%Ohio University EECS
%11-9-2010
%%  the code...
H = jacobian(sym_equations,sym_variables);
X = initial_estimate;

n_equations = 0;
if length(sym_equations)==length(sym_variables),
    n_equations = 1;
end

stop = 0;
while ~stop,
        F_X = subs(sym_equations,sym_variables,X);
        F_prime_X = subs(H,sym_variables,X);
        if ~isnumeric(F_prime_X),
            F_prime_X = eval(F_prime_X);
        end
    if n_equations ==1,
        d_X = (F_prime_X^-1)*F_X;
    else %overdetermined solution, use generalized inverse matrix
        d_X = ((F_prime_X.'*F_prime_X)^-1)*F_prime_X.'*F_X;
    end
    X = X - d_X.' ;
    if (sqrt(sum(d_X.^2)) < tolerance_rss),
        stop = 1;
    end
end
end
