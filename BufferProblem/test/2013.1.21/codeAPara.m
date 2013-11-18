%% Simulation Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.05;     % CFL condition
r_time      = 1/1000;  % Relaxation time
tEnd        = 0.1;      % End time
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 2;        % for NC = 1 , GH = 2
method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 1;        % IC: {1}Sod's, {2}LE, {3}RE, {4}DS, {5}SS, {6}Cavitation
plot_figs   = 1;        % 0: no, 1: yes please!
write_ans   = 1;        % 0: no, 1: yes please!
RRR=8.314;
% Using DG
P_deg       = 0;        % Polinomial Degree
Pp          = P_deg+1;  % Polinomials Points
% Using RK integration time step
RK_stages   = 4;        % Number of RK stages

gamma = 2.5; %1.4;        % Ratio of specific heats
 etpfix   = 0.90;	% {#} Harten's sonic entropy fix value, {0} no entropy fix
[MB_r,u,MB_p] = Euler_IC1d(x,1);
%% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    nv = 60;          % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;            % quadrature constant.
    w = w.*exp(v.^2); % weighting function of the Gauss-Hermite quadrature
    
    otherwise
        error('Order must be between 1 and 2');
end



%% Velocity-Space Grid:
% The actual nv value will be computed using 'lenght' vector function:
nv = length(v); 
% Using D.O.M.
    v = repmat(v,1,nx);     w = repmat(w,1,nx);
	p = zeros(nv,nx);   
    %% Load Macroscopic Velocity, Temperature and Fugacity
    [rho_0,u0,t0] = BGK_IC1d(x,IC_case);
 rho_0 = repmat(rho_0,nv,1);                       ux = repmat(u0,nv,1);       t = repmat(t0,nv,1);
    f0 = f_equilibrium_1d(rho_0,ux,v,t,RRR);

    %% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.
dt = dx*CFL/max(v(:,1)); 
dtdx = dt/dx;  % precomputed to save someflops

% Time domain discretization
time = 0:dt:tEnd;