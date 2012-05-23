%% 1D Heat Transfer FDM
% by Manuel Diaz
clc; clear; %close all;

%% Domain Grid
L = 30;    % Length
n = 30;    % Nodes
dx = L/n; dt = 0.5;
x  = 1:dx:L;

%% Initial Condition
fo = zeros(1,n);    %old
f  = zeros(1,n); 

%% Constants
alpha = 0.25;
iter = 400;         % time steps

%% Boundary Conditions
twall = 1;
f(1)  = twall;      %Dirichlet BC
fo(1) = twall;      %Dirichlet BC
f(n)  = f(n-1);     %Neumann BC
fo(n) = fo(n-1);    %Neumann BC

%% Main Loop
for k = 2:iter;
    for i = 2:n-1;
        f(i) = fo(i) + dt*alpha*(fo(i+1)-2*fo(i)+fo(i-1))/(dx^2);
    end
    fo = f;         %update domain
    fo(n) = f(n-1); %update BC
end

%% Make pretty figures
figure
plot(x,fo)
title '1D Heat Equation using FDM'
xlabel 'x cell'
ylabel 'Temperature'