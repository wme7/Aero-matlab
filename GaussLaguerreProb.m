%% Gauss Laguerre Quadrature
% solution for class exercise:

clear
clc
close all

%% Aproximate solution
% Tau(m) is a function of m.

m=5.5   % Exponent value of f(x).
n=9;    % number of points used to compute approximate solution.

[x w]=GaussLaguerre(n,0); % built in function to generate weight and points

% Initalizing row vectors:
l=length(x);
f=zeros(1,l);
t=zeros(1,l);

for i=1:l
    f(i)=x(i)^(m-1);
    t(i)=w(i)*f(i);
end

Gamma=sum(t)

%% Exact Solution
% Gamma=(m-1)! 
% Carefull! m can only be an integer number!

Gamma2=factorial(floor(m)-1) 
% I'm using a round to the floor in case m is a rational number.

%% Exact solution of the Gamma function
% this time using the Matlab's gamma function.

Gamma3=gamma(m)
