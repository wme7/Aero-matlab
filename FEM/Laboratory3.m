% Laboratory No.3
% Coded by Manuel Diaz, 2013.03.20
clear all; close all; clc;

%% Choose polynomial
pol = 2;
switch pol
    case{1} % polynomial #1
        P = @(x) 20 + x + 2*x.^2 + 3*x.^3 + 4*x.^4 + 5*x.^5 + 6*x.^6 + 7*x.^7;
    case{2} % polynomial #2
        P = @(x) 1./(1+x.^2);
end

%% define interval of integration
a = 1; 
b = 2*sqrt(3);

%% perform gauss integration
ngp = 1:5; % = {1,2,3,4,5}
approx = zeros(1,ngp);
for i = ngp
    approx(i) = gaussint(P,a,b,i);
end

%% compute exact solution
syms x;
exact = subs(int(P,x,a,b));

%% compute and plot error 
error = abs((exact-approx)/exact)*100;
plot(ngp,error,'or'); 
title('Integration error using n-gauss quadrature points')
ylabel('Error %'); xlabel('Number of gauss points')
