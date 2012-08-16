%% 1D Heat Transfer FDM
% by Manuel Diaz
clc; clear; close all;

%% Parameters
alpha = 0.25;       % difusion speed
iter = 400;         % iteration time steps

%% Domain Grid
L = 30;    % Length
n = 30;    % Nodes
dx = L/n; dt = 0.5;
x  = 1:dx:L;
t  = zeros(1,n);    % temperature of the nodes in our Domain

%% Wall Temperatura at x = 0
twall = 1;

%% Initial Condition
t_0 = zeros(1,n);   % @time = 0
t_0(1) = twall;     %Dirichlet BC
t_0(n) = t_0(n-1);   %Neumann BC

%% Main Loop
t_next = zeros(1,n); %Next time step
t = t_0;			 %Load I.C.
for k = 1:iter;
    for i = 2:n-1;
        t_next(i) = t(i) + (dt/dx^2)*alpha*(t(i+1)-2*t(i)+t(i-1));
    end
    % BC
    t_next(1) = twall;
	t_next(n) = t_next(n-1); 
	% Update info
	t = t_next;
end

%% Make pretty figures
plot(x,t,'.'); xlabel 'x cell'; ylabel 'Temperature'; title '1D Heat Equation using FDM';
