%**************************************************************************
%* Roe's Approximate 1D Riemann Solver
%**************************************************************************
%* Following the ideas of:
%* E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
%* Manchester U.K., Springer Editorial, 2nd Ed. Chapert 11.
%*
%* This code solves the Sod's shock tube problem 
%*
%* - t=0                               - t=tf
%* Density                             Density
%*   ****************|                 *********\
%*                   |                           \
%*                   |                            \
%*                   |                             ****|
%*                   |                                 |
%*                   |                                 ****|
%*                   ***************                       ***********
%*
%* coded by Manuel Diaz, 2012.12.25
%**************************************************************************
clear all; clc; % close all;
global gamma

%% Parameters
%cfl     = 0.9;     % Courant number
nx       = 100;     % number of cells
mx       = nx+1;    % number of nodes
ICx      = 1;      % IC: {1}~{12}. See Euler_IC1d.m
%tEnd    = 0.01;    % final time to compute
etpfix   = 0.90;	% {#} Harten's sonic entropy fix value, {0} no entropy fix
plot_fig = 1;       % {1} plot figures, {0} do NOT plot figures
wrt_sol  = 1;       % {1} write solution, {0} do NOT write solution file

%% Physical Constanst
%gamma = 1.4;         % Ratio of specific heats
gamma = 2;

%% Domain
xn = linspace(0,1,mx);          % cells nodes
dx = max(xn(2:mx)-xn(1:mx-1));  % Cell size
x  = xn(1:mx-1)+dx/2;        	% Cell centers

%% IC for Euler Riemann solver
[r,u,p,tEnd,cfl] = Euler_IC1d(x,ICx);

% Specific internal energy of the IC
e_s = p./((gamma-1).*r);
% Internal energy of the IC : e = e_s*rho
e = p./(gamma-1);
% Total Enegy: E = e + 1/2 rho*u^2
E = e + 1/2*r.*u.^2;
% Total Enthalpy: H = g/(g-1)p/rho + 1/2 u^2
H = gamma/(gamma-1)*p./r + u.^2/2;

%% Exact Riemann Solution
[xx,rhoexact,uexact,pexact,machexact,entroexact,energexact] = ...
   Exact_Riemann(r(1),u(1),p(1),r(nx),u(nx),p(nx),tEnd);

% %% Left and right conditions for Dirichlet BC's
%     r_left  = r(1);  % r_left
%     r_right = r(nx); % r_right
%     u_left  = u(1);  % u_left
%     u_right = u(nx); % u_right
%     p_left  = p(1);  % p_left
%     p_right = p(nx); % p_right
%     E_left  = E(1);  % E_left
%     E_right = E(nx); % E_right

%% Compute primitive Variables @ cell center variables
% Build vector U and dU, and F
    % Define:
    L = 1:nx-1; R = 2:nx; % inner nodes / middle points
    % U = [U1 U2 U3]' % defined at every cell center
    U = [r; r.*u; E];
        % U1 = Density, U2 = Momentum, U3 = Total Energy

% Initial time and Time step
t = 0;
        
%% Main Loop
n = 1;  % counter
U_next = zeros(3,nx);
%for time = t; 
while t < tEnd
    % Compute Roe Averages
    % Velovity 'u'
    u_bar = (sqrt (r(L)).*u(L) + sqrt (r(R)).*u(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Total Entalpy 'H'
    H_bar = (sqrt (r(L)).*H(L) + sqrt (r(R)).*H(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Sound Speed 'a'
    a_bar = sqrt((gamma-1)*(H_bar-0.5*u_bar.^2));

    % Compute Delta U's
    % dU = U(R)-U(L) at the cell boundaries { x_{i},x_{i+1}, ... }
    dU = U(:,R)-U(:,L);

    % Compute Fluxes at the cell centers
    % F = [F1 F2 F3]
    F = [r.*u; r.*u.^2 + p; r.*u.*H]; 
    % Fluxes to the left and right of 'F_ {i+1/2}'
    FL = F(:,L); FR = F(:,R);
    
    % Scaled Right Eigenvectors are given by:
    k1_bar = [1*ones(1,nx-1); u_bar-a_bar; H_bar-u_bar.*a_bar];
    k2_bar = [1*ones(1,nx-1); u_bar      ; 1/2*u_bar.^2      ];
    k3_bar = [1*ones(1,nx-1); u_bar+a_bar; H_bar+u_bar.*a_bar];
    
    % compute Roe waves strength alpha_bar
    alpha2_bar = (gamma-1)./(a_bar.^2).*(dU (1,:).*(H_bar-u_bar.^2) ...
        + u_bar.*dU(2,:) - dU(3,:));
    alpha1_bar = 1./(2*a_bar).*(dU (1,:).*(u_bar+a_bar) ...
        - dU(2,:)-a_bar.*alpha2_bar);
    alpha3_bar = dU(1,:)-(alpha1_bar + alpha2_bar);
    
    % Eigenvalues of A_bar are (same as the original A mat)
    lambda_bar(1,:) = abs(u_bar - a_bar);
    lambda_bar(2,:) = abs(u_bar);
    lambda_bar(3,:) = abs(u_bar + a_bar);
    
    % Entropy Fix
    boolv1 = lambda_bar(1,:) < etpfix;
        lambda_bar(1,:) = 0.5*(etpfix + lambda_bar (1,:).^2/etpfix).*boolv1 ...
            + lambda_bar(1,:).*(1-boolv1);
    boolv2 = lambda_bar(3,:) < etpfix;
        lambda_bar(3,:) = 0.5*(etpfix + lambda_bar (3,:).^2/etpfix).*boolv2 ...
            + lambda_bar(3,:).*(1-boolv2);
        
    % Conditioning data
    alpha1_bar = repmat(alpha1_bar,3,1);
    alpha2_bar = repmat(alpha2_bar,3,1);
    alpha3_bar = repmat(alpha3_bar,3,1);
    lambda1_bar = repmat(lambda_bar(1,:),3,1);
    lambda2_bar = repmat(lambda_bar(2,:),3,1);
    lambda3_bar = repmat(lambda_bar(3,:),3,1);

    % Update time step
    dt = dx*cfl/max(abs(lambda_bar(3,:)));
    t = t + dt;
    dtdx = dt/dx;
    
    % Roe Fluxes
    Flux = 0.5*(FL+FR)-0.5*(alpha1_bar.*lambda1_bar.*k1_bar + ...
                            alpha2_bar.*lambda2_bar.*k2_bar + ...
                            alpha3_bar.*lambda3_bar.*k3_bar );
                        
    % Compute next time step
    for i = 2:nx-1
        U_next(:,i) = U(:,i) - dtdx*(Flux(:,i) - Flux(:,i-1));
    end
    
    % BCs
    U_next(:,1)  = U(:,2); % Neumann condition to the left
    U_next(:,nx) = U(:,nx-1); % Neumann condition to the right
    
    % Compute variables of the new time step
    r_next = U_next(1,:);               % Density
    u_next = U_next (2,:)./U_next(1,:);  % Velocity
    E_next = U_next(3,:);               % Total Energy
    p_next = (gamma-1).*(E_next-r_next.*u_next.^2/2);  % Pressure
    e_next = 1/(gamma-1)*(p_next./r_next);      % Internal Energy
    a_bar_next = sqrt(gamma*p_next./r_next);    % sound speed
    m_next = u_next./a_bar_next;        % Mach 
    s_next = log(p_next./r_next.^gamma);% Entropy
    H_next = (E_next + p_next)./r_next; % Enthalpy
    
    % Update info
    U = U_next;
    r = r_next;
    u = u_next;
    e = e_next;
    p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
        
    % Update Counter
    n = n+1;
    
    % Plot figure
    if plot_fig == 1;
        subplot(2,3,1); plot(x,r,'o'); title('Velocity');
        subplot(2,3,2); plot(x,u,'o'); title('ernal Energy');
        subplot(2,3,3); plot(x,p,'o'); title('Pressure');
        subplot(2,3,4); plot(x,m,'o'); title('Mach number');
        subplot(2,3,5); plot(x,s,'o'); title('Entropy');
        subplot(2,3,6); plot(x,e,'o'); title('IntDensity');
    end
    drawnow
end

%% Write Results
if wrt_sol == 1;
    subplot(2,3,1); hold on; plot(x,r,'o'); title('Density');
    subplot(2,3,2); hold on; plot(x,u,'o'); title('Velocity');
    subplot(2,3,3); hold on; plot(x,p,'o'); title('Pressure');
    subplot(2,3,4); hold on; plot(x,m,'o'); title('Mach number');
    subplot(2,3,5); hold on; plot(x,s,'o'); title('Entropy');
    subplot(2,3,6); hold on; plot(x,e,'o'); title('Internal Energy');
    
    subplot(2,3,1); plot(xx,rhoexact);   hold off;
    subplot(2,3,2); plot(xx,uexact);     hold off;
    subplot(2,3,3); plot(xx,pexact);     hold off;
    subplot(2,3,4); plot(xx,machexact);  hold off;
    subplot(2,3,5); plot(xx,entroexact); hold off;
    subplot(2,3,6); plot(xx,energexact); hold off;
end
