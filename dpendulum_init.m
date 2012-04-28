%% Double Pendulum
% Clear work space and initaliza figure
clear; clc; close all; %figure;

%% Initial and Constant Values

% Initial values for tuned state (1):
% gamma0    = pi/12;
% dtgamma0  = pi;
% alpha0    = 0;
% dtalpha0  = 0;

% Initial values for tuned state (2):
% gamma0    = 0;
% dtgamma0  = 0;
% alpha0    = pi/6;
% dtalpha0  = pi/36;

% Initial values for tuned state (3):
% gamma0    = 0;
% dtgamma0  = 0;
% alpha0    = pi/24;
% dtalpha0  = pi/2;

% Initial values for chaotic solution:
% gamma0    = 0;
% dtgamma0  = 0;
% alpha0    = pi/2;
% dtalpha0  = 5; 

% Initial values for chaotic solution2:
gamma0    = 0;
dtgamma0  = 0;
alpha0    = pi/2;
dtalpha0  = 5.025;

% Constant Values
beta0     = pi/2; 
lambda    = 2.74; %rad/s
omega     = 5.48; %rad/s
psi       = 0.96; 
eta       = 0.24; 

fps       = 10;
%movie     = true;

%% Using solver ode45

duration  = 100;
ivp       = [gamma0; dtgamma0; alpha0; dtalpha0; beta0; lambda; omega; psi; eta];
[t,x]     = ode45(@dpendulum,[0 duration], ivp);

%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%[t,x]=ode45(@dpendulum,[0 duration], ivp, options);

%% Plot

nframes=duration*fps;
%t = linspace(0,duration,nframes);
%x = deval(sol,t);

gamma    = x(:,1);
dtgamma  = x(:,2);
alpha    = x(:,3);
dtalpha  = x(:,4);
% we assume:
l1=1; % b distance
l2=2; % c distance

%% Plot alpha vs Gamma
figure;
plot (t,x(:,3),'-.',t,x(:,1),'-')

%% Plot d_alpha vs d_Gamma
figure;
plot (t,x(:,4),'-.',t,x(:,2),'-')
