function xNew = Newton_n_system(f,df,x0)
%% Newton Method for N-system of equations
% Subroutine used to compute the newest for solving a nonlinear system of
% equations using Newton Method. Coded by Manuel Diaz, NTU, 2013.03.29.

% INPUT: f: function, df: system Jacobian, x0: Initial Guess.
% OUTPUT: xNew; new approximation.

delta_x = -df(x0)\f(x0); % solve for the increment
xNew = x0 + delta_x;     % add on to get new guess
%f(xNew)                 % uncoment to see if f(x) is really zero