%% One-dimensional Conservation Law Solver 
% Computing Burgers solution using a DG-FEM routine with Polynomial
% truncation for the nonlinear terms. 
%
% for solving the following problem:
%
% u_t + f(u)_x = s(u)
%
% where f(u) and S(u) can be
% For linear advection eq.: f(u) = a*u and s(u) = any function of u
% For non-linear advection: f(u) = u^2/2 and s(u) = any function of u 
%
% Function residual will be defined as:
%
% Residue(u) = - f(u)_x + s(u)
%
% Based on ideas of the following papers:
%
% 1. TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite Element
% Method for conservation laws II: General Framework. (1989)
% 2. Runge-Kutta Discontinuous Galerkin Method Using WENO Limiters. (2005)
%
% Coded by Manuel Diaz 2012.12.05

%% Clear Work Space 
clear all; close all; %clc;

%% Simulation Parameters
k         = 3;      % Space order / Number of degress of freedom: 0 to k
np        = k+1;    % Number of points per Cell/Element
quadn     = 3;      % element grid: {1}sLeg, {2}Lobatto, {3}Leg, {4}Radau
rks       = 3;      % Time order of the Sheme/RK stages
flux_type = 4;      % {1}Roe, {2}Global LF, {3}LLF, {4}Upwind (non-conservative)
a         = 1.0;    % scalar advection speed
cfl       = 1/(2*k+1);    % Courant Number
tEnd      = pi/15;  % Final Time for computation
nx        = 10;     % Number of Cells/Elements
MM        = 0.01;   % TVB constant M
IC_case   = 3;      % 4 cases are available
plot_figs = 1;      % {1}Plot figures, {0}Do NOT plot figures

%% Define Grid Cell's (Global) nodes
% Building nodes for cells/elements: 
x_left = 0; x_right = 1; dx = (x_right-x_left)/nx;
x_nodes = x_left : dx : x_right; % cells nodes

%% flux function
F = @(w) a * w;     % for Linear advection
%F = @(w) w.^2/2;    % for Invicid burger's eq.
% and derivate of the flux function
dF = @(w) a;        % for Linear advection
%dF = @(w) w;        % for Invicid burger's eq.

%% Source term fucntion
S = @(w) w.^2;       % for Source Term

%% SETUP
% 1. Build Cells/Elements (Local) inner points (quadrature points).
% 2. Build Weigthing values for our local.
% 3. Build Vandermonde Matrix for our local quadrature points.
[x,xi,w,V] = setup(k,x_nodes,quadn);

% Compute Math Objetcs:
if quadn == 1; bmath = 1; else bmath = 2; end;
switch bmath
    case{1} % Build Math objects for scaled Legendre polynomials. See Ref.[1]
        % M matrix
        M = diag([1 1/12 1/180 1/2800 1/44100 1/698544 1/11099088 1/176679360]);        
        % invM matrix
        invM = inv(M);
        % D matrix 
        D = [ ...
            0, 1, 0, (1/10), 0, (1/126), 0, (1/1716); ...
            0, 0, (1/6), 0, (1/70), 0, (1/924), 0; ...
            0, 0, 0, (1/60), 0, (1/756), 0, (1/10296); ...
            0, 0, 0, 0, (1/700), 0, (1/9240), 0; ...
            0, 0, 0, 0, 0, (1/8820), 0, (1/120120); ...
            0, 0, 0, 0, 0, 0, (1/116424), 0; ...
            0, 0, 0, 0, 0, 0, 0, (1/1585584); ...
            0, 0, 0, 0, 0, 0, 0, 0];
        % Scaled Legendre polynomials of deg 'l' evaluated at x = +1/2
        Ln = zeros(k+1,1); % column vector
        for l = 0:k; % Polynomials degree
        Ln(l+1) = sLegendreP(l,0.5);
        end
        % Scaled Legendre polynomials of deg 'l' evaluated at x = -1/2
        Lp = zeros(k+1,1); % column vector
        for l = 0:k; % Polynomials degree
        Lp(l+1) = sLegendreP(l,-0.5);
        end
        
    case{2} % Build Math objects for non-scaled Legendre polynomials See Ref.[3]
    l = (0:k)'; % Polynomials degree
    % M matrix
    M = diag(dx./(2*l+1));
    % invM matrix
    invM = inv(M);
    % D matrix 
    D = zeros(np,np);
    for ll = 0:k            % For all degress of freedom
        i = ll+1;           % Dummy index
        for j=1:np          % For all local points
            if j>i && rem(j-i,2)==1
                D(i,j)=2;   % D or diffentiated Legendre Matrix
            end
        end
    end
    % Scaled Legendre polynomials of degree 'l' evaluated at x = +1
    Ln = (1).^(l);  % LegP @ x_{i+1/2}^(-)
    % Scaled Legendre polynomials of degree 'l' evaluated at x = -1
    Lp = (-1).^(l); % LegP @ x_{i-1/2}^(+)
end

%% Load Initial Condition, u(x,0) = u0
u0 = u_zero(x,IC_case);
%f0 = F(u0);
%s0 = S(u0);

% Plot IC
% if plot_figs == 1; 
%     subplot(1,3,1); plot(x,u0); title('u0'); axis tight;
%     subplot(1,3,2); plot(x,f0); title('f0'); axis tight;
%     subplot(1,3,3); plot(x,s0); title('s0'); axis tight; 
% end;

%% Computing the Evolution of the residue 'L', du/dt = L(u) 
% Load Initial conditions
u = u0; %f = f0; s = s0;

% transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
% ut = zeros(np,nx);
% for l = 0:k             % for all degress of freedom
%     i = l+1;            % Dummy index
%     for j = 1:nx
%         ut(i,j) = (2*l+1)/2*sum(w(:,j).*u(:,j).*V(:,i));
%     end
% end
%     
% % transform f(x,t) to degress of freedom f(t)_{l,i} for each i-Cell/Element
% ft = zeros(np,nx); 
% for l = 0:k             % for all degress of freedom
%     i = l+1;            % Dummy index
%     for j = 1:nx
%         ft(i,j) = (2*l+1)/2*sum(w(:,j).*f(:,j).*V(:,i));
%     end
% end
% 
% % transform s(x,t) to degress of freedom s(t)_{l,i} for each i-Cell/Element
% st = zeros(np,nx);
% for l = 0:k             % for all degress of freedom
%     i = l+1;            % Dummy index
%     for j = 1:nx
%         %st(i,j) = (2*l+1)/2.*sum(w(:,j).*s(:,j).*V(:,i));
%         st(i,j) = dx/2*sum(w(:,j).*s(:,j).*V(:,i));
%     end
% end
% 
% % Volume term (test on methods)
% % 1. Cell Integral volume by Gauss-type quadrature.
% v_term1 = zeros(np,nx); 
% for l = 0:k
%     i = l+1;    % Dummy index
%     for j = 1:nx
%         v_term1(i,j) = (2*l+1)/2*sum(w(:,j).*f(:,j).*dsLegendreP(l,xi));
%     end
% end
% % 2. Cell Integral volume by T.W method
% v_term2 = D'*ft;

% Set Initial time step
t = 0; % time
n = 0; % counter
residue = zeros(np,nx);
%while t < tEnd
    % Time step 'dt'
    u_reshaped = reshape(u,1,nx*np); 
    dt  = dx*cfl/max(abs(u_reshaped));
    t  = t + dt;   % iteration time / increment time
    n  = n + 1;    % update counter
    
    % Plot solution every time step
    if plot_figs == 1; plot(x,u,'-'); axis([0,1,-1.5,1.5]); end;
    
    % Evaluate/update u, f and s terms on domain
    ut = V\u;
    f = F(u); ft = V\f;
    s = S(u); st = V\s;
    
    % Compute fluxes
    un = (ut'*Ln)'; % u_{i+1/2}^(-) -> Right u
    up = (ut'*Lp)'; % u_{i-1/2}^(+) -> Left u
    ub = [up(2:nx);un(1:nx-1)];
    h  = DGflux1d(F,dF,ub,flux_type); % Evaluate fluxes
     
    % BCs Periodic BC assuming element 1 and nx are ghost cells.
    ut(:,1)  = ut(:,nx-1);
    ut(:,nx) = ut(:,2);
    
    % Compute residue function:
    for i = 2:nx-1
        residue(:,i) = ( D'*ft(:,i) - ...          % Volume term
                        h(i)*Ln + h(i-1)*Lp + ...  % flux terms
                        st(:,i) ).*diag(M);        % Source term
    end
     
    % Compute next time step
    ut_next = ut + dt*residue;

    % UPDATE info
    ut = ut_next;
    
    % Update plot
    drawnow

    % Transform degress u(t)_{l,i} into values u(x,t)
    u = (ut'*V')'; 
    
%end % time loop 

%% Examine Output
if plot_figs == 1; 
u = (ut'*V')'; % f = F(u); % s = S(u);
figure
    subplot(2,1,1); plot(x,u,'o-'); title('u-final'); axis tight;
    subplot(2,1,2); plot(x,residue,'o-'); title('residue'); axis tight;
end;