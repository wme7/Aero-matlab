function [f_s,f_r] = Upwindflux1d(u,a)
% Traditional Upwinding Method for scalar advection equation.
% coded by Manuel Diaz, NTU, 2012.07.20

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
f_s = zeros(1,nx); %fsl = zeros(1,nx); fsh = zeros(1,nx); 
f_r = zeros(1,nx); %frl = zeros(1,nx); frh = zeros(1,nx); 
    
for j = x
    % Compute fluxes for TVD
    f_r(j) = a_p(j)*u(j) + a_m(j)*u(j+1);
    f_s(j) = a_p(j)*u(j-1) + a_m(j)*u(j);
end

return