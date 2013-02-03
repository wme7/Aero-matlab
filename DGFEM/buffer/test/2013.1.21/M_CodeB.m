function [rhor,uxr,tr,pr,fr] = ModEuler(h,f0,rho0,ux0,t0,v)
%% Global Variables
global CFL r_time theta dt dtdx nx
global w k nv 
global gamma etpfix

%% Split ICs in R and L
% fr and fl
    fr = (1-h).*f0;         fl = h.*f0;
% rho_r and rho_l
    rhor = (1-h).*rho0;     rhol = h.*rho0;
% u_r and u_l
    ur = (1-h).*ux0;        ul = h.*ux0;
% t_r and t_l
    tr = (1-h).*t0;         tl = h.*t0;

%% Compute new equilibriums for R and L sides
% theta = 0;
Mfr_eq = f_equilibrium_1d(rhor,ur,v,tr,theta);  % right
Mfl_eq = f_equilibrium_1d(rhol,ul,v,tl,theta);  % left
Mf_eq  = f_equilibrium_1d(rho0,ux0,v,t0,theta); % total

%% Physical Constanst
gamma = 2.5; %1.4;  % Ratio of specific heats
etpfix   = 0.90;	% {#} Harten's sonic entropy fix value, {0} no entropy fix

e_s = MB_p./((gamma-1).*MB_r);
% Internal energy of the IC : e = e_s*rho
e = MB_p./(gamma-1);
% Total Enegy: E = e + 1/2 rho*u^2
E = e + 1/2*MB_r.*u.^2;
% Total Enthalpy: H = g/(g-1)MB_p/rho + 1/2 u^2
H = gamma/(gamma-1)*MB_p./MB_r + u.^2/2;

 %% Left and right conditions for Dirichlet BC's
  
    r_left  = MB_r(1);  % r_left
    r_right = MB_r(nx); % r_right
    u_left  = u(1);  % u_left
    u_right = u(nx); % u_right
    p_left  = MB_p(1);  % p_left
    p_right = MB_p(nx); % p_right
    E_left  = E(1);  % E_left
    E_right = E(nx); % E_right

%% Compute primitive Variables @ cell center variables
% Build vector U and dU, and F
    % Define:
    L = 1:nx-1; R = 2:nx; % inner nodes / middle points
    % U = [U1 U2 U3]' % defined at every cell center
    MB_r=(1-h).*MB_r ;  E=(1-h).*E;
    U = [MB_r; MB_r.*u; E];
    
        % U1 = Density, U2 = Momentum, U3 = Total Energy
        
        % Compute Roe Averages
    % Velovity 'u'
    u_bar = (sqrt (MB_r(L)).*u(L) + sqrt (MB_r(R)).*u(R)) ...
                ./ (sqrt(MB_r(L))+sqrt(MB_r(R)));
    % Total Entalpy 'H'
    H_bar = (sqrt (MB_r(L)).*H(L) + sqrt (MB_r(R)).*H(R)) ...
                ./ (sqrt(MB_r(L))+sqrt(MB_r(R)));
    % Sound Speed 'a'
    a_bar = sqrt((gamma-1)*(H_bar-0.5*u_bar.^2));

    % Compute Delta U's
    % dU = U(R)-U(L) at the cell boundaries { x_{i},x_{i+1}, ... }
    dU = U(:,R)-U(:,L);

    % Compute Fluxes at the cell centers
    % F = [F1 F2 F3]
    F = [MB_r.*u; MB_r.*u.^2 + MB_p; MB_r.*u.*H]; 
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
    
    % Roe Fluxes
    Flux = 0.5*(FL+FR)-0.5*(alpha1_bar.*lambda1_bar.*k1_bar + ...
                            alpha2_bar.*lambda2_bar.*k2_bar + ...
                            alpha3_bar.*lambda3_bar.*k3_bar );
                       [M1,M2,M3] =M_ROE(k,w,f,v,h);
                       MMM=[M1;M2;M3];
                       dMMM = MMM(:,R)-MMM(:,L);
                        for i =2 : nx-1                            
                            U_next(:,i) = U(:,i) - dtdx*(1-h(1,i))*(Flux(:,i) - Flux(:,i-1))-dtdx*dMMM(:,i);
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
    MB_r = r_next;
    u = u_next;
    e = e_next;
    MB_p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
        
    % Update Counter
    n = n+1;
        
        
        
        
        
        
        
        
        
        
        
        