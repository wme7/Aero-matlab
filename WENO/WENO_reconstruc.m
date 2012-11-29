% WENO 1D reconstruction for equation 
%
% u_t + f_x = 0 (advection equation in conservation form)
% 
% Where f(u) = a*u(x), assuming "a" is on not constant along domain x.

% Compute size of domain:
nx = length(u);

%% Test varialble "a"  
% Check whether a and b are scalar velocities or 
% prescribed velocities in the entire domain
[m,n,r] = size(a);
if m == 1 && n == 1 && r == 1
    a = a*ones(1,nx); % map the x-velocity 
else
    %do nothing
end
clear m n r;

%% Flux Spliting
% Using Upwinding flux spliting:
up = 1/2*(a.*u + max(a).*u); % flux{+}
un = 1/2*(a.*u - max(a).*u); % flux{-}

%% compute Fluxes in the x Domain using:
% WENO3

% WENO5

%% Compute Flux: u_i+1/2
u_flux = uup + uun;