%% RATE OF CONVERGENCE OF LINEAR FINITE ELEMENT METHOD FOR HELMHOLTZ EQUATION
%
% This example is to show the rate of convergence of linear finite element
% approximation of the Helmholtz equation on the unit square:
%
% $$- \Delta u + k^2 u = f \; \hbox{in } (0,1)^2.$$
%
% for the following boundary condition:
%

%% 
clear all; close all;
global k        % wave number
[node,elem] = squaremesh([0,1,0,1],0.5);
%[node,elem]  =  PMLdomain(node,h);
% pde = Helmholtzdata2; %sommerfeld radiation condition

%% Options
option.L0 = 2;
option.printlevel = 1;
option.maxIt = 4;
option.plotflag = 1;
option.solver = 'shifted';
option.shifted = 1;

%% Dirichlet boundary condition
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Helmholtzdata1;
wavenumberlist = [1 5 10 40 100];
for i = 1:length(wavenumberlist)
    k = wavenumberlist(i); 
    femHelmholtz(node,elem,pde,bdFlag,option);
    display(k);
    %%
end
%% 
% As the wavenumber increases, the discretization is less stable. For small
% k, optimal order of convergence is still observed but for larger k, one
% needs small enough h to recovery the optimal order of convergence. Note
% that the exact solution is smooth, the loss of accuracy is from the
% stability of discretization.

%% Absorbing boundary condition
bdFlag = setboundary(node,elem,'ABC');
pde = Helmholtzdata2;
wavenumberlist = [1 5 10 40 100];
for i = 1:length(wavenumberlist)
    k = wavenumberlist(i); 
    femHelmholtz(node,elem,pde,bdFlag,option);
    display(k);
    %%
end
%%
% We change the solution to be sin(k*pi*x)^2*sin(k*pi*y)^2 with frequency
% depending on k. Then the deterior of rate is more obvious since now the
% approximability requires a small h.