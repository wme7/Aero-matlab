function [f_s,f_r] = TVDflux1d(u,a,dtdx)
%% Grid Size
nx=length(u);
x = 2:nx-1;

%% Test a  
% Check whether a and b are scalar velocities or 
% prescribed velocities in the entire domain
[m,n,r] = size(a);
if m == 1 && n == 1 && r == 1
    a = a*ones(1,nx); % map the x-velocity 
else
    %do nothing
end
clear m n r;

% Compute Flux Spliting Factors
a_p = max(a,0); % a{+}
a_m = min(a,0); % a{-}

%% Compute TVD fluxes
% Note:
% s = Left flux & r = Right flux
% l = low flux  & h = high flux

% Initalize Arrays
fsl = zeros(1,nx); fsh = zeros(1,nx); f_s = zeros(1,nx);
frl = zeros(1,nx); frh = zeros(1,nx); f_r = zeros(1,nx);
    
for j = x
    % Compute fluxes for TVD
    frl(j) = a_p(j)*u(j) + a_m(j)*u(j+1);
%     frh(j) = (1/2)*a(j)*(u(j)+u(j+1)) - (1/2)*(a(j)^2)*dtdx*(u(j+1)-u(j));
%     f_r(j) = frl(j) + phi_x(j)*( frh(j) - frl(j) );
    f_r(j) = frl(j);
    
    
    fsl(j) = a_p(j)*u(j-1) + a_m(j)*u(j);
%     fsh(j) = (1/2)*a(j)*(u(j-1)+u(j)) - (1/2)*(a(j)^2)*dtdx*(u(j)-u(j-1));
%     f_s(j)  = fsl(j) + phi_x(j-1)*( fsh(j) - fsl(j) );
    f_s(j)  = fsl(j)
end

return