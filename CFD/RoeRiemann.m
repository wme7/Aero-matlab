%**************************************************************************
%* The Riemann solver of Roe 
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
clear all; close all; clc

%% Parameters
cfl      = 0.1;     % courant number
nx       = 40;      % number of cells
n        = nx-1;    % number of nodes at boundary
tEnd     = 0.05;    % final time to compute
plot_fig = 1;       % {1} plot figures, {2} do NOT plot figures

%% Physical Constanst
gamma = 1.4; 
Cv = 1.4182;
CP = gamma*Cv;
R_gas = 8.314;

%% Domain
dx = (1-0)/nx;
x = 0+dx/2:dx:1-dx/2;

%% IC for Euler Riemann solver
% Sod's Problem
%x = [xL , xR  ];
r0 = [1 , 0.125]; 
u0 = [0 , 0    ];
p0 = [1 , 0.100];

% Pre-Allocate variables
r = zeros(1,nx); u = zeros(1,nx); p = zeros(1,nx);

for i = 1:nx %for all cells 
    if x(i) < 0.5
        r(i) = r0(1); 
        u(i) = u0(1);
        p(i) = p0(1);
    else
        r(i) = r0(2);
        u(i) = u0(2);
        p(i) = p0(2);
    end
end

% Exact speed of sound
a = sqrt(gamma*p./r);

% System's exact eigenvalues
lambda(1,:) = u - a;
lambda(2,:) = u;
lambda(3,:) = u + a;

% Use exact eigenvalues to compute initial time step 'dt'.
lambda_bar = lambda;

t = 0;  % time 
m = 0;  % counter
U_next = zeros(3,nx);
while t < tEnd
%for i = 1:20
    dt = dx*cfl/max(max(abs(lambda_bar)));
    t = t + dt; % time
    m = m + 1;  % counter
    clear lambda_bar;
    
    % Plot figure
    if plot_fig == 1;
        subplot(2,2,1); plot(x,r,'.'); axis tight; title('density');
        subplot(2,2,2); plot(x,u,'.'); axis tight; title('velocity');
        subplot(2,2,3); plot(x,p,'.'); axis tight; title('pressure');
    end
        
    % Energy
    e = p./((gamma-1)*r);
    % Total Energy
    E = e + u.^2/2;
    % Total Entalpy
    H = (E + p)./r;
    
    % Build vector U and dU, and F
    % Define:
    L = 1:nx-1; R = 2:nx; % inner nodes
    % U = [u1 u2 u3]' % defined at every cell
    U = [r; r.*u; r.*E];
    % dU = U(R)-U(L) at the cell boundaries ( x_{i+1/2} )
    dU = U(:,R)-U(:,L);
    % F = [u1 u2 u3]' % defined at every cell
    F = [r.*u; r.*u.^2+p; r.*u.*(E+p./r)];
    
    % Compute Roe Averged values:
    % Velovity 'u'
    u_bar = (sqrt(r(L)).*u(L) + sqrt(r(R)).*u(R)) ...
        ./ (sqrt(r(L))+sqrt(r(R)));
    % Total Entalpy 'H'
    H_bar = (sqrt(r(L)).*H(L) + sqrt(r(R)).*H(R)) ...
        ./ (sqrt(r(L))+sqrt(r(R)));
    % Sound Speed 'a'
    a_bar = (gamma-1)*(H_bar-1/2.*u_bar.^2/2).^(1/2);
    
    % Eigenvalues of A_bar are (same as the orginal A mat)
    lambda_bar(1,:) = u_bar - a_bar;
    lambda_bar(2,:) = u_bar;
    lambda_bar(3,:) = u_bar + a_bar;
    
    % Scaled Right Eigenvectors are given by:
    k1_bar = [1*ones(1,nx-1); u_bar-a_bar; H_bar-u_bar.*a_bar];
    k2_bar = [1*ones(1,nx-1); u_bar      ; 1/2*u_bar.^2      ];
    k3_bar = [1*ones(1,nx-1); u_bar+a_bar; H_bar+u_bar.*a_bar];
    
    % compute Roe waves strength alpha_bar
    alpha2_bar = (gamma-1)./(a_bar.^2).*(dU(1,:).*(H_bar-(u_bar).^2) ...
        + u_bar.*dU(2,:) - dU(3,:));
    alpha1_bar = 1./(2*a_bar).*(dU(1,:).*(H_bar-a_bar) ...
        - dU(2,:)-a_bar.*alpha2_bar);
    alpha3_bar = dU(1,:)-(alpha1_bar + alpha2_bar);
    
    % compute Flux_{i+1/2}
    % conditioning data
    FL = F(:,L); FR = F(:,R);
    alpha1_bar = repmat(alpha1_bar,3,1);
    alpha2_bar = repmat(alpha2_bar,3,1);
    alpha3_bar = repmat(alpha3_bar,3,1);
    lambda1_bar = repmat(lambda_bar(1,:),3,1);
    lambda2_bar = repmat(lambda_bar(2,:),3,1);
    lambda3_bar = repmat(lambda_bar(3,:),3,1);
    
    % compute flux in cells boundaries
    Flux = 1/2*(FL+FR)-1/2*(alpha1_bar.*abs(lambda1_bar).*k1_bar + ...
                            alpha2_bar.*abs(lambda2_bar).*k2_bar + ...
                            alpha3_bar.*abs(lambda3_bar).*k3_bar );
    
    % Compute next time step
    for i = 2:nx-1
        U_next(:,i) = U(:,i) + dt/dx*(Flux(:,i)-Flux(:,i-1));
    end
    
    % BCs
    U_next(:,1) = U(:,2); % Neumann condition at the left
    U_next(:,nx) = U(:,nx-1); % Neumann condition at the right
    
    % Compute Variables of new time step
    r_next = U_next(1,:);               % Density
    u_next = U_next(2,:)./U_next(1,:);  % Velocity
    E_next = U_next(3,:)./U_next(1,:);  % Total Energy
    e_next = E_next - u_next.^2;        % specific internal energy
    p_next = (gamma-1).*r.*e_next;      % pressure
    
    % Update info
    r = r_next;
    u = u_next;
    p = p_next;
    
    % plot
    drawnow

end