% HEQDEMO  This program creates the H-equation example in Chapter 2.
% Solve the H-equation with the default parameters in nsold and plot
% the results.
% 
%
global A_heq;
n=100;
c=.9;
%
% Set the nodal points for the midpoint rule.
%
mu=1:n; mu=(mu-.5)/n; mu=mu';
%
% Form and store the kernel of the integral operator in a global variable.
%
cc=.5*c/n;
A_heq=ones(n,1)*mu'; A_heq=cc*A_heq'./(A_heq+A_heq');
%
tol=[1.d-6,1.d-6];
x=ones(n,1);
%
% Use the default parameters in nsold.m
%
[hc, errsd, ierrd]=nsold(x, 'heq', tol);
%
% Plot the results.
%
plot(mu,hc);
xlabel('\mu'); ylabel('H','Rotation',1);
