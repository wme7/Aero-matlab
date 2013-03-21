% Using our function gauss1d to integrate a given polynomials
% Coded by Manuel Diaz 2013.03.20
clear all; close all; clc;

%-------------------------------------------------------------------------
% (a.) define P(x):
P = @(x) 20 + x + 2*x.^2 + 3*x.^3 + 4*x.^4 + 5*x.^5 + 6*x.^6 + 7*x.^7;

% compute gauss locations and weights
[w,psi] = gauss1d(5);

% compute the change the integration range to -1 to 1
a = 1; b = 2*sqrt(3);
x = (b+a)/2 + (b-a)/2*psi;

% perform quadrature
results  = (b-a)/2*sum(w.*P(x));

% display my resutl
disp(results);

%-------------------------------------------------------------------------
%(b.) define Q(y)
Q = @(y) 1./(1+y.^2);

% compute gauss locations and weights
[W,Psi] = gauss1d(5);

% compute the change the integration range to -1 to 1
a = 1; b = 2*sqrt(3);
y = (b+a)/2 + (b-a)/2*Psi;

% perform quadrature
results2  = (b-a)/2*sum(W.*Q(y));

% display my 2nd result
disp(results2);