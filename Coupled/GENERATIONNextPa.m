% Compute Semiclassical ICs
% if (input ~= 7 && input ~= 8)
    
    E = p_total(1,:)+(0.5).*r_total(1,:).*u_total(1,:).^2; % Energy
    t = 4*E./r_total(1,:)-2*u_total(1,:).^2;    % Temperature
    r =r_total(1,:)./sqrt(pi*t);    % Fugacity
    
    
     r = repmat(r,nv,1); ux = repmat(u_total,nv,1); t = repmat(t,nv,1);
    
% else
%     % Do nothing
% end

%% Load Selected case Initial condition:
% number of points required
 nx = length(x);

% Preallocating
r_0 = zeros(1,nx); 
u_0 = zeros(1,nx); 
t_0 = zeros(1,nx); 

% Parameters of regions dimensions
x_middle = ceil(nx/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;

% Initial Condition for our 2D domain
% Fugacity
r_0(l_1) = r(1); % region 1
r_0(l_2) = r(2); % region 2
% Velovity in x
u_0(l_1) = u(1); % region 1
u_0(l_2) = u(2); % region 2
% temperature
t_0(l_1) = t(1); % region 1
t_0(l_2) = t(2); % region 2