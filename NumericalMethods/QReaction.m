%% Reaction Solver
% Homework 7, Problem 2
clc
close all
clear all
%% Convencion
% c(1) = Cs; c(2) = Ce, c(3) = Ces, c(4) = Cp

%% Initial Conditions
c = [1 5.0e-5 0.0 0.0];

%% Time Range
ti=0
tf=1e2

%% Solving with ode45
options = odeset('AbsTol',1e-6);
[T,C] = ode45(@reaction_func,[ti tf],c,options);

% Make figure
figure
loglog(T,C(:,1),'blue',T,C(:,2),'red',T,C(:,3),'green',T,C(:,4),'black');
axis([0 1e2 1e-10 1e1]);
legend('Cs','Ce','Ces','Cp');
set(gca,'fontsize',12);
grid on

%% Solving with ode23s
[T2,C2] = ode23s(@reaction_func,[ti tf],c);
% Make figure
figure
loglog(T2,C2(:,1),'blue',T2,C2(:,2),'red',T2,C2(:,3),'green',T2,C2(:,4),'black');
axis([0 1e2 1e-10 1e1]);
legend('Cs','Ce','Ces','Cp');
set(gca,'fontsize',12);
grid on

%% Discusion
%
% By comparing the computer time used by solving these ode system, we
% observer that there is a substantial difference when using ode45 and 
% ode23s to solve a stiff system of equations.
%
% Such time is evidenced by observing the grafical output results and
% comparing the number of points used computed to calculate a solution for
% the same time range:
%
%   - ode45 need 236253 points
%   - ode23s need 22 points