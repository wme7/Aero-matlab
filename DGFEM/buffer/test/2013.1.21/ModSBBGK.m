function [fr_next] = ModSBBGK(h,h_next,M,fr,Ml,v,flxtype)
%function [rhor,uxr,pr,Er,tr,zr,fr_next] = ModSBBGK(h,h_next,M,fr,Ml,v,flxtype)
%% Global Variables
global CFL r_time theta dt dtdx nx %#ok<NUSED>
global w k nv 
global gamma etpfix %#ok<NUSED>

% M_eq = M;
% v = c;

%% Main Loop
% initialize variables
fr_next = zeros(1,nx);
u_eq = zeros(1,nx);
u_a  = zeros(1,nx);
u_b  = zeros(1,nx);

for i = 1:nv
    % load subcase
    u_eq(:) = M(i,:);
    u_a(:) = Ml(i,:);
    u_b(:) = fr(i,:);
    
    % Compute TVD Fluxes
    [Fa_left,Fa_right] = Upwindflux1d(u_a,v(i,:));
    [Fb_left,Fb_right] = Upwindflux1d(u_b,v(i,:));
    
    % Compute next time step
    u_b_next = u_b - dtdx*h(i,:).*(Fa_right - Fa_left) ...
                        - dtdx*h(i,:).*(Fb_right - Fb_left) ...
                        +(dt/r_time).*h(i,:).*(u_eq - (fr(i,:) + Ml(i,:))) ...
                        +(fr(i,:) + Ml(i,:)).*(h(i,:) - h_next(i,:));
    % BC
    u_b_next(1) = u_b_next(2);
    u_b_next(nx) = u_b_next(nx-1);
       
    % Going back to f
    fr_next(i,:) = u_b_next(:); % time step n+1 or n+1/2
end