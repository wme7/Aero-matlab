function [tout, yout] = euler(rhs, t0, y0, h, N)

% euler('rhs',t0,y0,h,N)
%
% integrates a single ordinary differential equation
% of the form 
%     y'=f(t,y)
% where 
%   the initial value of t = t0
%   and y(t0) = y0 
%   the stepsize = h
%   and the number of steps = N
% Note: the right hand side must be defined in the file rhs.m

% allocate space and initialize to vectors to zero
tout = zeros(N+1,1);
yout = zeros(N+1,1);

% set initial values
tout(1) = t0;
yout(1) = y0;

% apply Euler's method
for i=1:N
  tout(i+1) = tout(i) + h;
  slope = feval(rhs,tout(i),yout(i));
  yout(i+1) = yout(i) + slope*h;
end
