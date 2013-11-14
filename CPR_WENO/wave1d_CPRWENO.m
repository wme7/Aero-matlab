%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with CPR/FR
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2013.11.13
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% 1. A flux reconstruction approach to high-order schemes including
% Discontinuous Galerkin methods. by H.T. Huynh, AIAA 2007.
% 2. A simple essentially non-oscillatory limiter for the correction
% procedure vie reconstruction (CPR) framework. by Du, Shu & Zhang 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 
% 1. Basic CPR scheme implementation without RK integration method.
% 2. Implementation of the 'simple 'WENO limiter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Parameters
fluxfun = 'nonlinear'; % select flux function
cfl = 0.02; % CFL condition
tEnd = 1; % final time
K = 3; % degree of accuaracy
nE = 20; % number of elements
M = 100; % MODminmod parameter

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=1; flux = @(w) a*w; 
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
    
    %% Limiting
    
    % Build u_j(x) lagrange polynomials
    ulp = zeros(size(x,2),size(x,1));
    for j = 1:nE
        ulp(j,:) = lagrange(x(:,j),u(:,j)); % interpolation
    end
    
    % Build Cell averages for every E_j
    u_bar = w*u/2;
    
    % detect troubled cells
    tcd = TroubleCellDectector(u_bar,'MODminmod',1,u(1,:),u(1+K,:),M,dx);
    tCells = tcd.troubledCells; disp(tCells);
    
    term0 = zeros(1,size(ulp,2));
    term1 = zeros(1,size(ulp,2));
    term2 = zeros(1,size(ulp,2));
    for j = tCells
        if ( j~=1 && j~=nE )
            % Build smooth indicators
            for s = 1:K
                dpl0 = ulp(j-1,:);
                dpl1 = ulp(j,:);
                dpl2 = ulp(j+1,:);
                % Derivate 's' times
                for i = 1:s
                    dpl0 = polyder(dpl0);
                    dpl1 = polyder(dpl1);
                    dpl2 = polyder(dpl2);
                end
                % Do: dx^(2s-1)*integrate(dp^2,a,b)
                integ0 = dx^(2*s-1)*polyint(conv(dpl0,dpl0));
                term0(s) = polyval(integ0,x(1+K,j-1))-polyval(integ0,x(1,j-1));
                integ1 = dx^(2*s-1)*polyint(conv(dpl1,dpl1));
                term1(s) = polyval(integ1,x(1+K,j))-polyval(integ1,x(1,j));
                integ2 = dx^(2*s-1)*polyint(conv(dpl2,dpl2));
                term2(s) = polyval(integ2,x(1+K,j+1))-polyval(integ2,x(1,j+1));
            end
            Beta0 = sum(term0);
            Beta1 = sum(term1);
            Beta2 = sum(term2);
            
            % Build Weights
            gamma = [1e-6,0.999998,1e-6]; epsilon = 1e-6;
            w_tilde0 = gamma(1)./(epsilon+Beta0).^2;
            w_tilde1 = gamma(2)./(epsilon+Beta1).^2;
            w_tilde2 = gamma(3)./(epsilon+Beta2).^2;
            wsum = w_tilde0 + w_tilde1 + w_tilde2;
            w0 = w_tilde0./wsum;
            w1 = w_tilde1./wsum;
            w2 = w_tilde2./wsum;
            
            % Modify troubled polynomials
            P0 = ulp(j-1,:);
            P1 = ulp( j ,:);
            P2 = ulp(j+1,:);
            P0_bar = zeros(size(P0));
            P1_bar = zeros(size(P1));
            P2_bar = zeros(size(P2));
            P0_bar(K+1) = u_bar(j-1);
            P1_bar(K+1) = u_bar(j);
            P2_bar(K+1) = u_bar(j+1);
            P0_tilde = P0 - P0_bar + P1_bar;
            P2_tilde = P2 - P2_bar + P1_bar;
            u_new = w0*P0_tilde + w1*P1 + w2*P2_tilde;
            u(:,j) = polyval(u_new,x(:,j)); %update info
        end
    end
            
    %% CPR
    
    % compute fluxes at node coordinates
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
    %u_nface(1) = u_nface(end); % left BD
    %u_pface(end) = u_pface(1); % right BD
    
    % Apply Neumann BCs
    u_nface(1) = u_pface(1); % left BD
    u_pface(end) = u_nface(end); % right BD
    
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
    
    % Plot u
    plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
        
    %if rem(it,10) == 0
        drawnow;
    %end
    
end