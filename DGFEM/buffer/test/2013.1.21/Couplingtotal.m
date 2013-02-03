%% Coupling Case 
% Shck Tube Problem
% Kinetic equation use BGK and Hydrodynamic equation use Roe Euler
% Yun-da Tsai   
% Reference Manuel Diaz
clear all; clc; % close all;

%%  Parameters
name        ='SBBGK1d';        % Simulation Name
CFL         = 0.05;            % CFL condition
tEnd        = 0.1;             % End time
method      = 1;               % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 1;               % IC: {1}Sod's, {2}LE, {3}RE, {4}DS, {5}SS, {6}Cavitation
plot_figs   = 1;               % 0: no, 1: yes please!
write_ans   = 1;               % 0: no, 1: yes please!
RK_stages   = 4;               % Number of RK stages
%%  Parameters-BGK  Part
r_time      = 1/10000;         % Relaxation time
theta       = 0;               % {-1} BE, {0} MB, {1} FD.
quad        = 2;               % for NC = 1 , GH = 2
%% Space Discretization
nx       = 100;                 % number of cells
mx       = nx+1;                % number of nodes
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x
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
%%  Parameters-Roe  Euler Part
etpfix   = 0.90;	            % {#} Harten's sonic entropy fix value, {0} no entropy fix

dt       = dx*CFL/max(v(:,1));
dtdx     = dt/dx; 
gamma    = 2.5;                 % Ratio of specific heats

%% Define a ID name for results file
% [ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,r_time,IC_case);
%% Open a Files to store the Results
% if write_ans == 1
%     file = fopen(IDn,'w');
%     % 'file' gets the handel for the file "case.plt".
%     % 'w' specifies that it will be written.
%     % similarly 'r' is for reading and 'a' for appending.
%     fprintf(file, 'TITLE = "%s"\n',ID);
%     fprintf(file, 'VARIABLES = "x" "density" "velocity" "energy" "pressure" "temperature" "fugacity"\n');
% end

%% Velocity-Space Grid:
% The actual nv value will be computed using 'lenght' vector function:
nv = length(v); 

% Using D.O.M.
    v = repmat(v,1,nx);     w = repmat(w,1,nx);

% Initialize Arrays
%    ux = zeros(nv,nx);      t = zeros(nv,nx);       
%    r = zeros(nv,nx);       n = zeros(nv,nx);
	p = zeros(nv,nx);   
    %% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
    [r0,u0,t0] = SSBGK_IC1d(x,IC_case);
    
% Using Discrete Ordinate Method:
    r = repmat(r0,nv,1); ux = repmat(u0,nv,1); t = repmat(t0,nv,1);

% Compute distribution IC of our mesoscopic method by assuming the equilibrium 
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
    f0 = f_equilibrium_1d(r,ux,v,t,theta);

% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1
   figure(1)
   surf(f0); grid on;
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end
% Compute Initial Macroscopic Momemts:
    [n,j_x,E] = macromoments1d(k,w,f0,v);
%% IC for Euler Riemann solver
ICx=1;
[r,u,p]=Euler_IC1d(x,ICx);


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
%% Left and right conditions for Dirichlet BC's
r_left  = r(1);  % r_left
r_right = r(nx); % r_right
u_left  = u(1);  % u_left
u_right = u(nx); % u_right
p_left  = p(1);  % p_left
p_right = p(nx); % p_right
E_left  = E(1);  % E_left
E_right = E(nx); % E_right
 %% Compute primitive Variables @ cell center variables
% Build vector U and dU, and F
% Define:
L = 1:nx-1; R = 2:nx; % inner nodes / middle points
% U = [U1 U2 U3]' % defined at every cell center
U = [r; r.*u; E];
% U1 = Density, U2 = Momentum, U3 = Total Energy
%% Cut Function
for i = 2:nx-1
if (x(1,i) < 0.335)
     h(1,i)=1;
   elseif( x(1,i) > 0.675)
     h(1,i)=0;
   else
   h(1,i)= (x(1,i)-0.675)/(0.335-0.675);
   end
end
%% Main Loop 

time = 0:dt:tEnd;
% Compute next time step
for tsteps = time 
    for i = 2:nx-1

        
        
        if (x(1,i) < 0.335)
         % compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(r,ux,v,t,theta);
            
            % initialize variables
            u_next_BGK = zeros(1,nx);
            u_eq_BGK = zeros(1,nx);
            u_BGK = zeros(1,nx);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                u_eq_BGK(:) = f_eq(i,:);
                u_BGK(:) = f(i,:);
                
                % Compute the smoothness factors, r(j), from data, u(j).
                [r] = theta1d(u_BGK,a(i));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(u_BGK,a(i),dtdx,phi);

                % Compute next time step
                u_next_BGK = u_BGK - dtdx*(F_right - F_left) ...
                    + (dt/r_time)*(u_eq_BGK-u_BGK);

                % BC
                u_next_BGK(1) = u_next_BGK(2);
                u_next_BGK(nx) = u_next_BGK(nx-1);

                % UPDATE info
                u_BGK = u_next_BGK;
                
                % Going back to f
                f_BGK(i,:) = u_BGK(:);
            end
            
            % Compute macroscopic moments
            [n,j_x,E] = macromoments1d(k,w,f_BGK,v);
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing chalenge)
            [r,ux,t,p] = macroproperties1d(n,j_x,E,nx,nv,theta);
            
           
    
  a(1,i) = n;
  b(1,i) = j_x;
  c(1,i) = E  ;
  

 
        elseif   ( x(1,i) > 0.675)
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
    
    % Roe Fluxes
    Flux = 0.5*(FL+FR)-0.5*(alpha1_bar.*lambda1_bar.*k1_bar + ...
                            alpha2_bar.*lambda2_bar.*k2_bar + ...
                            alpha3_bar.*lambda3_bar.*k3_bar );

    U_next_R(:,i) = U_R(:,i) - dtdx*(1-h(1,i))*(Flux(:,i) - Flux(:,i-1))
  
    % BCs
    U_next_R(:,1)  = U_R(:,2);     % Neumann condition to the left
    U_next_R(:,nx) = U_R(:,nx-1); % Neumann condition to the right
    
    % Compute variables of the new time step
    r_next = U_next_R(1,:);               % Density
    u_next = U_next_R (2,:)./U_next_R(1,:);  % Velocity
    E_next = U_next_R(3,:);                       % Total Energy
    p_next = (gamma-1).*(E_next-r_next.*u_next.^2/2);  % Pressure
    e_next = 1/(gamma-1)*(p_next./r_next);      % Internal Energy
    a_bar_next = sqrt(gamma*p_next./r_next);    % sound speed
    m_next = u_next./a_bar_next;                % Mach 
    s_next = log(p_next./r_next.^gamma);        % Entropy
    H_next = (E_next + p_next)./r_next;         % Enthalpy
    
    % Update info
    U = U_next_R;
    r = r_next;
    u = u_next;
    e = e_next;
    p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
        
    % Update Counter
    n = n+1;
    
  

else
           
            % compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(r,ux,v,t,theta);
            
            % initialize variables
            u_next_C = zeros(1,nx);
            u_eq_C = zeros(1,nx);
            u_C = zeros(1,nx);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                u_eq_C(:) = f_eq(i,:);
                u(:) = f_C(i,:);
                
                % Compute the smoothness factors, r(j), from data, u(j).
                [r] = theta1d(u,a(i));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(u,a(i),dtdx,phi);

                % Compute next time step
                u_next_C = u_C - dtdx*(F_right - F_left)-dtdx*...
                    (u_eq_C-u_eq_C) + (dt/r_time)*(u_eq_C-u_C);

                % BC
                u_next_C(1) = u_next_C(2);
                u_next_C(nx) = u_next_C(nx-1);

                % UPDATE info
                u_C = u_next_C;
                
                % Going back to f
                f_C(i,:) = u_C(:);
            end
            
            % Compute macroscopic moments
            [n,j_x,E] = macromoments1d(k,w,f,v);
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing chalenge)
            [r,ux,t,p] = macroproperties1d(n,j_x,E,nx,nv,theta);
            
           d(1,i) = n;
           e(1,i) = j_x;
           f(1,i) = E  ;
%    ------------------------------------------------------------------------        
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
    F_C = [r.*u; r.*u.^2 + p; r.*u.*H]; 
    % Fluxes to the left and right of 'F_ {i+1/2}'
    FL_C = F(:,L); FR_C = F(:,R);
    
    % Scaled Right Eigenvectors are given by:
    k1_bar_C = [1*ones(1,nx-1); u_bar-a_bar; H_bar-u_bar.*a_bar];
    k2_bar_C = [1*ones(1,nx-1); u_bar      ; 1/2*u_bar.^2      ];
    k3_bar_C = [1*ones(1,nx-1); u_bar+a_bar; H_bar+u_bar.*a_bar];
    
    % compute Roe waves strength alpha_bar
    alpha2_bar_C = (gamma-1)./(a_bar.^2).*(dU (1,:).*(H_bar-u_bar.^2) ...
        + u_bar.*dU(2,:) - dU(3,:));
    alpha1_bar_C = 1./(2*a_bar).*(dU (1,:).*(u_bar+a_bar) ...
        - dU(2,:)-a_bar.*alpha2_bar);
    alpha3_bar_C = dU(1,:)-(alpha1_bar + alpha2_bar);
    
    % Eigenvalues of A_bar are (same as the original A mat)
    lambda_bar_C(1,:) = abs(u_bar - a_bar);
    lambda_bar_C(2,:) = abs(u_bar);
    lambda_bar_C(3,:) = abs(u_bar + a_bar);
    
    % Entropy Fix
    boolv1 = lambda_bar_C(1,:) < etpfix;
        lambda_bar_C(1,:) = 0.5*(etpfix + lambda_bar_C (1,:).^2/etpfix).*boolv1 ...
            + lambda_bar_C(1,:).*(1-boolv1);
    boolv2 = lambda_bar_C(3,:) < etpfix;
        lambda_bar_C(3,:) = 0.5*(etpfix + lambda_bar_C (3,:).^2/etpfix).*boolv2 ...
            + lambda_bar_C(3,:).*(1-boolv2);
        
    % Conditioning data
    alpha1_bar_C = repmat(alpha1_bar_C,3,1);
    alpha2_bar_C = repmat(alpha2_bar_C,3,1);
    alpha3_bar_C = repmat(alpha3_bar_C,3,1);
    lambda1_bar_C = repmat(lambda_bar_C(1,:),3,1);
    lambda2_bar_C = repmat(lambda_bar_C(2,:),3,1);
    lambda3_bar_C = repmat(lambda_bar_C(3,:),3,1);
    
    % Roe Fluxes
    Flux_C = 0.5*(FL_C+FR_C)-0.5*(alpha1_bar_C.*lambda1_bar_C.*k1_bar_C + ...
                            alpha2_bar_C.*lambda2_bar_C.*k2_bar_C + ...
                            alpha3_bar_C.*lambda3_bar_C.*k3_bar_C );

    U_next_C(:,i) = U_C(:,i) - dtdx*(1-h(1,i))*(Flux_C(:,i) - Flux_C(:,i-1))
  
    % BCs
    U_next_C(:,1)  = U_C(:,2);     % Neumann condition to the left
    U_next_C(:,nx) = U_C(:,nx-1); % Neumann condition to the right
    
    % Compute variables of the new time step
    r_next = U_next_C(1,:);               % Density
    u_next = U_next_C (2,:)./U_next_C(1,:);  % Velocity
    E_next = U_next_C(3,:);                       % Total Energy
    p_next = (gamma-1).*(E_next-r_next.*u_next.^2/2);  % Pressure
    e_next = 1/(gamma-1)*(p_next./r_next);      % Internal Energy
    a_bar_next = sqrt(gamma*p_next./r_next);    % sound speed
    m_next = u_next./a_bar_next;                % Mach 
    s_next = log(p_next./r_next.^gamma);        % Entropy
    H_next = (E_next + p_next)./r_next;         % Enthalpy
    
    % Update info
    U = U_next_C;
    r = r_next;
    u = u_next;
    e = e_next;
    p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
   
 
 U_next_C + [n j_x E]
   end
end
end
        
    
    






