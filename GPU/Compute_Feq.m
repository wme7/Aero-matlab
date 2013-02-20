%% Evaluating Equilibrium Distribution Function
clear all; close all; clc;

% Parameters
IC_case = 5;
D = gpuDevice(1);
quad = 2;
theta = 0;
GPU = 1;

%% Space Discretization
nx = 400; ny = 400; % for 6Gb of memory 435 by 435 is the limit
X = linspace(0,1,nx); dx = abs(max(X(2:nx)-X(1:nx-1)));
Y = linspace(0,1,ny); dy = abs(max(Y(2:ny)-Y(1:ny-1)));
[x,y] = meshgrid(X,Y);

%% Velocity Discretization:
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-10,10];  % range: a to b
    nv = 50;        % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % cotes Degree 5

    case{2} % Gauss Hermite Quadrature:
    nv = 20;        % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
     k = 1;         % quadrature constant.
     w = w.*exp(v.^2);
    
    otherwise
        error('Order must be between 1 and 2');
end
tic
%% Velocity-Space Grid:
% The actual nv value will be computed using 'length' vector function:
nv = length(v); 
if GPU == 1
     v = gpuArray(v); w = gpuArray(w);
end

% Using Discrete Ordinate Method (Replicating Data)
    [vx,vy] = meshgrid(v,v');       [wx,wy] = meshgrid(w,w');
    vx = repmat(vx,[1,1,ny,nx]);    wx = repmat(wx,[1,1,ny,nx]);
    vy = repmat(vy,[1,1,ny,nx]);    wy = repmat(wy,[1,1,ny,nx]);


%% Initial Conditions
% Semiclassical ICs: Fugacity[z], Velocity[u] and Temperature[T] 
    [z0,ux0,uy0,t0] = SSBGK_IC2d(x,y,IC_case);

    % load to GPU memory
if GPU ==1
     z0 = gpuArray(z0);  ux0 = gpuArray(ux0);
     uy0 = gpuArray(uy0);  t0 = gpuArray(t0);
     theta = gpuArray(theta);
end

% Using Discrete Ordinate Method: (Replicating Data for easyer arrays Ops)
z = reshape(z0,1,1,ny,nx);    z = repmat(z,[nv,nv,1]);
ux = reshape(ux0,1,1,ny,nx);  ux = repmat(ux,[nv,nv,1]);
uy = reshape(uy0,1,1,ny,nx);  uy = repmat(uy,[nv,nv,1]);
t = reshape(t0,1,1,ny,nx);    t = repmat(t,[nv,nv,1]);

time(1) = toc

tic
% Compute distribution IC of our mesoscopic method by assuming the equilibrium
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
if GPU == 1
    f0_gpu = arrayfun(@f_equilibrium_2d,z,ux,uy,t,vx,vy,theta);
else
    f0 = f_equilibrium_2d(z,ux,uy,t,vx,vy,theta);
end
time(2) = toc

f = f0_gpu;

% gather f0 to host
tic
if GPU ==1
    f0 = gather(f0_gpu);
end
time(3) = toc

% Total Time
totaltime = sum(time)

%% Figure
% Plot Probability Distribution Function, f, for an expecific point in space:
figure(1)
i = 2; j = 2;
zz = f0(:,:,j,i);
surf(zz); grid on;
xlabel('v1 - Velocity Space in x');
ylabel('v2 - Velocity Space in y');
zlabel('f - Probability');
clear zz i j

%% Clear GPU Memory
% reset(D);


