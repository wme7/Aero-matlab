% Testing ENO Stencils
clear all; close all; clc;
 
% Compute coefcients of polynomial of order k=3
k = 3; 
for i = 1:k+1 % dummy index
    r = i-2;  % r = -1:2
    for l = 1:3  % dummy index
        j = l-1; % j = 0:2
        c(i,l) = c_rj(r,j,k);
    end
end

%% Create an IC
x = 0:0.1:1;
% u = 0.5 + sin(2*pi*x); % sine
u = 0*(x>0.5) + (x<0.5); % riemann

% Compute flux at cell points x_i
f = 2*u;

%% Compute fluxes at cell interfaces x_i+1/2
% convention: x_{i-r+j} where r is the left shift
i = 4; 
% r = -1, j = 0:2 
g1 = c(1,1)*u(i+1) + c(1,2)*u(i+2) + c(1,3)*u(i+3); fprintf('g1: %1.6f\n',g1);
% r = 0 , j = 0:2
g2 = c(2,1)*u( i ) + c(2,2)*u(i+1) + c(2,3)*u(i+2); fprintf('g2: %1.6f\n',g2);
% r = 1 , j = 0:2
g3 = c(3,1)*u(i-1) + c(3,2)*u( i ) + c(3,3)*u(i+1); fprintf('g3: %1.6f\n',g3);
% r = 2 , j = 0:2
g4 = c(4,1)*u(i-2) + c(4,2)*u(i-1) + c(4,3)*u( i ); fprintf('g4: %1.6f\n',g4); 
