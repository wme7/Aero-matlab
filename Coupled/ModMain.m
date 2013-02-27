%-------------------------------------------------------------------------
% Coupling Kinetic and Fluid System of Euler Equations
% To solve Shock Tube Problem
%
% Kinetic equation use SBBGK and Hydrodynamic equation use Roe Euler
%
% Based on: 
%  [1] Pierre Degond, Giacomo Dimarco and Luc Mieussens
%      A moving interface method for dynamic kinetic-fluid coupling
%      Journal of Computation Physics 227(2007)1176-1208
%
% By Manuel Diaz and Yun-Da Tsai   
% 007@IAM 25.01.2013
%-------------------------------------------------------------------------

clear all; clc; close all;

%% Global Variables
global CFL r_time theta dt dtdx nx
global w k nv
global gamma etpfix

%% Controling Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.05;     % CFL condition
r_time      = 1/10000;  % Relaxation time
tEnd        = 0.1;      % End time
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 2;        % for NC = 1 , GH = 2
method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 1;        % IC: {1}Sod's, {2}LE, {3}RE, {4}DS, {5}SS, {6}Cavitation
plot_figs   = 1;        % 0: no, 1: yes please!
write_ans   = 1;        % 0: no, 1: yes please!
gamma       = 3.8;      % Ratio of specific heats
flxtype     = 2;        % {1} Roe, {2} LF, {3} LLF, {4} Upwind <-non-conservative!
etpfix      = 0.90;     % {#} Harten's sonic entropy fix value, {0} no entropy fix

%% Space Discretization
nx  = 100;                      % number of cells
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Load Initial Condition
[z0,ux0,t0,p0,rho0,E0] = SSBGK_IC1d(x,IC_case);

%% Load Initial Cut Function
xa = 0;   % buffer left boundary
xb = 1;   % buffer right boundary
 h0  = cutfunc(x,xa,xb);  % Physical Cut Function
%  h0 = ones(size(x));
% h0 = zeros(size(x));

%% Discretization of the Velocity Space
% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [c,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    nv = 60;          % nodes desired (the actual value)
    [c,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;            % quadrature constant.
    w = w.*exp(c.^2); % weighting function of the Gauss-Hermite quadrature
    
    otherwise
        error('Order must be between 1 and 2');
end

%% Applying DOM
% The actual nv value will be computed using 'lenght' vector function:
nv = length(c); 
% Remap velocity points
    c = repmat(c,1,nx);     w = repmat(w,1,nx);
% Remap classical IC
    [rho0,ux0,p0] = apply_DOM(rho0,ux0,p0,nv);
% Remap Semiclassical IC
    [z0,t0,E0] = apply_DOM(z0,t0,E0,nv);
% Remap h coeficient
  h0 = repmat(h0,nv,1);

%% Semi0-classical Equilibrium Distribution Function
M0 = f_equilibrium_1d(z0,ux0,c,t0,theta); 

%% Load initial Conditions and Spliting of information
M = M0; rho = rho0; ux = ux0; t = t0; p = p0; z = z0; h = h0; E = E0;

%% Split ICs in R and L
rhor = h.*rho0;     rhol = (1-h).*rho0;
uxr = h.*ux0;       uxl = (1-h).*ux0;
pr = h.*p0;         pl = (1-h).*p0;

[zr,~,tr,~] = macroproperties1d(rhor,rhor.*uxr,pr+0.5*rhor.*uxr.^2,nx,nv,theta);
[zl,~,tl,~] = macroproperties1d(rhol,rhol.*uxl,pl+0.5*rhol.*uxl.^2,nx,nv,theta);

fr = f_equilibrium_1d(zr,uxr,c,tr,theta);
Ml = f_equilibrium_1d(zl,uxl,c,tl,theta) ;

Ml(isnan(Ml)) = 0;     fr(isnan(fr)) = 0;
f = fr + Ml; % computed here for ploting purposes

%% Main Loop 
% Compute next time step
dt = dx*CFL/max(c(:,1));
time  = 0:dt:tEnd;
dtdx = dt/dx;

count = 1; % iteration counter
tic;
for tsteps = time
      
    % Compute M_eq for the entire domain
    M_eq = f_equilibrium_1d(z,ux,c,t,theta);     % total % theta = 0;
    
    % Plot IC
    if plot_figs == 1; surf(f); end;
    
    % Compute vector 'q'
    q = [rhol(1,:) ; rhol(1,:).*uxl(1,:) ; pl(1,:)+0.5*rhol(1,:).*uxl(1,:).^2];
    
    % Update Physical cut function 'h'
    h_next = h; % this means: fixed buffer assumption

% uncomment for BGK ******************************************************
    %%  Evaluate Modified Boltzmann BGK
    [fr_next] = ModSBBGK(h,h_next,M_eq,fr,Ml,c,flxtype);
    
    % Compute macroscopic moments
%     [rhor,rhour,Er] = macromoments1d(k,w,fr_next,c);
    
    % macroscopic properties
%     [zr,uxr,tr,pr] = macroproperties1d(rhor,rhour,Er,nx,nv,theta);

%         rhol = zeros(size(x));
%         rhoul = zeros(size(x));  
%         El = zeros(size(x));
%         Ml_next = zeros(size(fr));
% uncomment for BGK ******************************************************
    
% uncomment for EULER*****************************************************   
%     rhor = zeros(size(x));
%     rhour = zeros(size(x));
%     Er = zeros(size(x));
%     fr_next = zeros(size(Ml));

    %%  Evaluate Modificed Roe Euler Solver
    [rhol,rhoul,El] = ModEuler(h,h_next,q,fr_next,c);
    
    %plot partial result
%     if plot_figs == 1; 
%         subplot(1,3,1); plot(x,rhol,'o'); title('Density');
%         subplot(1,3,2); plot(x,rhoul./rhol,'o'); title('Velocity');
%         subplot(1,3,3); plot(x,El,'o'); title('Energy');
%     end
%     
    % macroscopic properties
    [zl,uxl,tl,pl] = macroproperties1d(rhol,rhoul,El,nx,nv,theta);
    
    % Apply DOM
    [tl,zl,uxl] = apply_DOM(tl,zl,uxl,nv);
    
    % Compute Ml_next
    Ml_next = f_equilibrium_1d(zl,uxl,c,tl,theta) ;
    
%%  uncomment for EULER***************************************************** 
    
    % New time step info: sum left and righ values with 'NaN' filter sum
    % function. 
    Ml_next(isnan(Ml)) = 0;     fr_next(isnan(fr_next)) = 0;
    f_next   = fr_next + Ml_next;    % total f
    
    [rho_next,rhou_next,E_next] = macromoments1d(k,w,f_next,c);

    [z_next,ux_next,t_next,p_next] = macroproperties1d(rho_next,rhou_next,E_next,nx,nv,theta);
    [rho_next,ux_next,p_next] = apply_DOM(rho_next,ux_next,p_next,nv);
    [t_next,z_next,E_next]= apply_DOM(t_next,z_next,E_next,nv);
    
    
    
    
      rhol_next = (1-h).*rho_next;
       uxl_next = (1-h).*ux_next;
        pl_next = (1-h).*p_next;
        
        rhor_next = h.*rho_next;     
        uxr_next = h.*ux_next;       
        pr_next = h.*p_next;   
    
    [zl_next,~,tl_next,~] = macroproperties1d(rhol_next,rhol_next.*uxl_next,pl_next+0.5*rhol_next.*uxl_next.^2,nx,nv,theta);
    [zr_next,~,tr_next,~] = macroproperties1d(rhor_next,rhor_next.*uxr_next,pr_next+0.5*rhor_next.*uxr_next.^2,nx,nv,theta);
    
    
    fr_next = f_equilibrium_1d(zr_next,uxr_next,c,tr_next,theta);
    Ml_next = f_equilibrium_1d(zl_next,uxl_next,c,tl_next,theta) ;
     
 

    
    
    
    
    
%      rho = nansum([rhol;rhor]);  % total density
%      rhou= nansum([rhoul;rhour]);% total velocity in x
%      E   = nansum([El;Er]);      % total persure
%     
    %Apply Neumann BC's in total variables
    f_next(:,1)  = f_next(:,2);           f_next(:,end)  = f_next(:,end-1);
    rho_next(:,1)= rho_next(:,2);                   rho_next(:,end)= rho_next(:,end-1);

        % Recover Semiclassical Conditions for next time strhep
  
    
    % Update information
      f= f_next;
    fr = fr_next;
    rho=rho_next;
    p=p_next;
    ux=ux_next;
    z=z_next;
    t=t_next;
    Ml=Ml_next;
    
    rhol=rhol_next;
     uxl=uxl_next;
     p=pl_next;
    
    
    
    % update counter
    count = count+1;
    
    % update plot
    drawnow
end
toc;
% write/plot final output
%plot(x,r_total,'o');








