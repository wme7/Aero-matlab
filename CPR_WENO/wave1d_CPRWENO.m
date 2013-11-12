%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with CPR/FR
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: A flux reconstruction approach to high-order schemes including
% Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 1.Basic Scheme Implementation without RK integration method.
% 2. Basic implementation of the WENO limiter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.01; % CFL condition
tEnd = 1; % final time
K = 5; % degree of accuaracy
nE = 40; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=-1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' %Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 1],nE,'LGL',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; 
x = xgrid.nodeCoordinates; w = xgrid.weights';

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);

% Build Lagrange k-Polynomials
l = LagrangePolynomial(xgrid.solutionPoints);
L.lcoef = double(subs(l.lagrangePolynomial,-1));
L.rcoef = double(subs(l.lagrangePolynomial,1));
L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));

% IC
u0 = IC(x,3);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),0.9*min(min(u0)),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; it = 0;

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    % Plot u
    %plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
    
    % compute fluxes in node coordinates
    f = flux(u);

    % Interpolate u and flux values at the boundaries of Ij
    switch xgrid.quadratureType
        case 'LGL'
            u_lbd = u(1,:);
            u_rbd = u(end,:);
            f_lbd = f(1,:);
            f_rbd = f(end,:);
        otherwise
            u_lbd = L.lcoef*u;
            u_rbd = L.rcoef*u;
            f_lbd = L.lcoef*f;
            f_rbd = L.rcoef*f;
    end
    % Build Numerical fluxes acroos faces
    u_pface = [u_lbd,0]; % + side 
    u_nface = [0,u_rbd]; % - side 

    % Apply Periodic BCs
    u_nface(1) = u_nface(end); % left BD
    u_pface(end) = u_pface(1); % right BD

    % LF numerical flux
    alpha = max(max(abs(dflux(u)))); 
    nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
    nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

    % flux derivate
    df = L.dcoef*f;

    % Compute the derivate: F = f + gL*(nfluxL-f_bdL) + gR*(nfluxR-f_bdR)
    dF = df + dg.RR*(nfluxL - f_lbd) + dg.RL*(nfluxR - f_rbd);

    % next time info!
    u_next = u - dt*dF/J;
    
    % update info
    u = u_next;
    
    if rem(it,10) == 0
        drawnow;
    end
%end

%% Find Troubled cells:
% Build Cell averages for every E_j
%u0_bar = w*u0/dx;
u_bar = w*u/2; M = 100;

u_1tilde = u_nface(1:end-1) - u_bar;
u_2tilde = u_bar - u_pface(2:end);

% Taking into account periodic BCs!
Dpu_bar = [u_bar(2:end),u_bar(1)] - u_bar;
Dnu_bar = [u_bar(end),u_bar(1:end-1)] - u_bar;

% modif to take care of unsigned zeros
A = [u_1tilde;Dpu_bar;Dnu_bar]; A(find([u_1tilde;Dpu_bar;Dnu_bar]==0)) = 1e-16; 
B = [u_2tilde;Dpu_bar;Dnu_bar]; B(find([u_2tilde;Dpu_bar;Dnu_bar]==0)) = 1e-16;

uMOD_1tilde = MODminmod(A,M,dx);
uMOD_2tilde = MODminmod(B,M,dx);

% Mark troubled cells:
troubleE = find(uMOD_1tilde ==0 & uMOD_2tilde == 0);

%% WENO reconstruction for troubled cells
% Compute Smooth indicators
Bcoef = zeros(K,K+1);
for s = 1:K 
    syms x;
    dpdx  = l.dnlagrangePolynomial(s);
    Bcoef(s,:) = double(int((dx)^(2*s-1)*(dpdx).^2,x,-1,1));
end

% Beta factors for every element
B = sum(Bcoef*u);

%%
gamma = [1e-6,0.999998,1e-6];
epsilon = 1e-6;
w_tilde(1) = gamma(1)./(epsilon+B(troubleE)).^2;
w_tilde(2) = gamma(2)./(epsilon+B(troubleE)).^2;
w_tilde(3) = gamma(3)./(epsilon+B(troubleE)).^2;
w0 = w_tilde(1)/sum(w_tilde);
w1 = w_tilde(2)/sum(w_tilde);
w2 = w_tilde(3)/sum(w_tilde);
toc % 0.001 sec

switch K
    case 3
        phi0 = LGL_K3(xi+2);
        phi1 = LGL_K3(xi);
        phi2 = LGL_K3(xi-2);
    case 4
        phi0 = LGL_K4(xi+2);
        phi1 = LGL_K4(xi);
        phi2 = LGL_K4(xi-2);
    case 5
        phi0 = LGL_K5(xi+2);
        phi1 = LGL_K5(xi);
        phi2 = LGL_K5(xi-2);
end

u_new2 = zeros(size(u));
for j = troubleE
    P0 = phi0*u(:,j-1);
    P1 = phi1*u(:,j);
    P2 = phi2*u(:,j+1);
    P0_bar = u_bar(j-1);
    P1_bar = u_bar(j);
    P2_bar = u_bar(j+1);
    P0_tilde = P0 - P0_bar + P1_bar;
    P2_tilde = P2 - P2_bar + P1_bar;
    u_new2(:,j) = w0*P0_tilde + w1*P1 + w2*P2_tilde;
end

%update info
u = u_new2;

end
