%INITIAL_CONDITION
% Generates initial conditions for u,v and p
% Functions called: ubd, vbd

tic

if geval == 1		% Horizontal Poiseuille flow to the right
%    u0 = (YU - yv(1)).*(yv(end) - YU)*6/(yv(end)-yv(1))^3; % Average = 1
%   Integration of velocity profile with Simpson's rule
    u0 = ((YU - DYU/2 - yv(1)).*(yv(end) - YU + DYU/2) +...
         4*(YU - yv(1)).*(yv(end) - YU) +...
         (YU + DYU/2 - yv(1)).*(yv(end) - YU - DYU/2))...
         /(yv(end)-yv(1))^3; 
    v0 = zeros(size(XV));
    p0 = (xu(end) - XP)*12/(Re*(yv(end)-yv(1))^2);
elseif geval == 2		% Vertical Poiseuille flow
    u0 = zeros(size(XU));
%   Integration of velocity profile with Simpson's rule
    v0 = ((XV - DXV/2 - xu(1)).*(xu(end) - XV + DXV/2) +...
          4*(XV - xu(1)).*(xu(end) - XV ) +...
	  (XV + DXV/2 - xu(1)).*(xu(end) - XV - DXV/2))...
	  /(xu(end) -xu(1))^3; % Average = 1
    p0 = (yv(end) - YP)*12/(Re*(xu(end) -xu(1))^2);
elseif geval == 3		% Backward facing step
    u0 = zeros(size(XU)); v0 = zeros(size(XV)); p0 = zeros(size(XP));
    k = K + 1 - ny(end):K;
    u0(k,:) = 6*(YU(k,:) - yseglen(end)).*(yv(end) - YU(k,:))...
         /(yv(end) - yseglen(end))^3;
elseif geval == 4		% Driven cavity
    u0 = zeros(size(XU)); v0 = zeros(size(XV)); p0 = zeros(size(XP));
elseif geval == 5		% Uniform flow under angle alpha
    u0 = cos(alpha)*ones(size(XU)); v0 = sin(alpha)*ones(size(XV)); 
    p0 = ones(size(XP));
elseif geval == 6		% Uniform flow under angle alpha
    u0 = cos(alpha)*ones(size(XU)); v0 = sin(alpha)*ones(size(XV)); 
    p0 = ones(size(XP));
elseif geval == 7		% Horizontal Poiseuille flow to the left
%    u0 = -(YU - yv(1)).*(yv(end) - YU)*6/(yv(end)-yv(1))^3; % Average = 1
%   Integration of velocity profile with Simpson's rule
    u0 = -((YU - DYU/2 - yv(1)).*(yv(end) - YU + DYU/2) +...
         4*(YU - yv(1)).*(yv(end) - YU) +...
         (YU + DYU/2 - yv(1)).*(yv(end) - YU - DYU/2))...
         /(yv(end)-yv(1))^3; 
    v0 = zeros(size(XV));
    p0 = -(xu(end) - XP)*12/(Re*(yv(end)-yv(1))^2);
else
  error('Wrong value for geval in INITIAL CONDITION')
end
u0 = u0'; u0 = u0(:); v0 = v0'; v0 = v0(:); p0 = p0'; p0 = p0(:);

tijd = toc; disp(['initial_condition time = ',num2str(tijd)])
