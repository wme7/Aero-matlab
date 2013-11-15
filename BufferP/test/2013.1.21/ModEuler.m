function [rhol,rhoul,El] = ModEuler(h,h_next,q,fr,v)
%% Global Variables
global CFL r_time theta dt dtdx nx %#ok<*NUSED>
global w k nv 
global gamma etpfix

%% Compute primitive Variables @ cell center variables
% Build vector U and dU, and F
    % Define:
    L = 1:nx-1; R = 2:nx; % inner nodes / middle points
    % here U = [U1 U2 U3]' is defined for every cell center
    U = q;
        % U1 = Density, U2 = Momentum, U3 = Total Energy
    
    r = q(1,:);
    u = q(2,:)./q(1,:);
    E = q(3,:);
    p = (gamma-1).*(E-r.*u.^2/2);
    H = (E + p)./r;

%% Main Loop
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

%% Compute next time step
% compute difference of the integral(v*fl)
[M1,M2,M3] =integrate_vf(k,w,fr,v); M = [M1;M2;M3]; dM = M(:,R) - M(:,L);

% Remap from h(nx,nx) to h(3,nx)
hh = h(1:3,:);
hh_next = h_next(1:3,:);
dhh = hh_next - hh;

% Compute next time step
U_next = zeros(3,nx);
for i = 2 : nx-1
    U_next(:,i) = U(:,i) - dtdx.*(1-hh(:,i)).*(Flux(:,i) - Flux(:,i-1)) - ...
                           dtdx.*(1-hh(:,i)).* dM(i)- ...
                           q(:,i).*dhh(:,i);
end

% BCs
 U_next(:,1)  = U(:,2); % Neumann condition to the left
 U_next(:,nx) = U(:,nx-1); % Neumann condition to the right

% Compute variables of the new time step
rhol = U_next(1,:);     % Density
rhoul = U_next (2,:);   % Momentum
El = U_next(3,:);       % Total Energy
