% Integration using generalized Gauss-Laguerre
% coded by Manuel Diaz, 2012.04.01
clear all; close all; clc;

% paramaters
nu = 1/2;
z = 0.35;

% fucntion
x = 0:0.05:15;
f = @(x)(exp(x).*x.^nu)./(exp(x)-z).^2;

% Using gauss-laguerre quadrature
tic;
[xi,wi] = GaussLaguerre(20,nu);
wi = wi.*exp(xi).*xi.^-nu;
q_GL = sum(wi.*f(xi));
toc

% Integrating with matlab
tic;
q_exact = quad(f,0,15,1e-12);
toc

% plot
plot(x,f(x),'.')

% Error
disp('Max Error:');
disp(q_GL-q_exact);
