function [F] = calcLagrangianForce(Xl, L, dtheta, XT, K, Kp);
%
% [F] = calcLagrangianForce(Xl, L, dtheta, XT, K, Kp);
%
%  Calculates the elastic force at each Lagrangian Point
%  
%  Returns:
%     F  = elastic force on each Lagrangian point
%    
%  Input:
%     Xl     = Lagrangian point positions
%     L      = number of Lagrangian points
%     dtheta = Lagrangian point spacing
%     XT     = target point config
%
%  Notes:
%     Does not handle wrapping across periodic boundaries of
%     the immersed boundary. 
%
%
%  License: This code is free to use for any purposes, provided
%           any publications resulting from the use of this code
%           reference the original code/author.
%
%  Author:  Samuel Isaacson (isaacson@math.utah.edu)
%  Date:    11/2007
%
%  Please notify the author of any bugs, and contribute any
%  modifications or bug fixes back to the original author.
%
%  Disclaimer:
%   This code is provided as is. The author takes no responsibility 
%   for its results or effects.



persistent Lap isInit L0Up L0Dn idxUp idxDn tmp LUp LDn TUp TDn numL;


if( size( isInit ) == 0 )
  isInit = 1; 
  
  e = ones(L,1);
  Lap = spdiags([e -2*e e], [-1:1], L, L);
  
  % make periodic:
  Lap(1,L) = 1;
  Lap(L,1) = 1;
  
  Lap = (K / (dtheta * dtheta)) * Lap;
  
  numL = length(XT(:,1));
  
  idxUp = [2:numL 1]';
  idxDn = [numL 1:(numL-1)]';
  tmp     = XT(idxUp,:) - XT;
  L0Up    = sqrt(dot(tmp,tmp,2));
  tmp     = XT - XT(idxDn,:);
  L0Dn    = sqrt(dot(tmp,tmp,2));
  
  LUp     = zeros(numL,1);
  LDn     = zeros(numL,1);
  TUp     = zeros(numL,2);
  TDn     = zeros(numL,2);
end

% no rest length springs
F = [( Lap * Xl(:,1) ) ( Lap * Xl(:,2) )]; 


% springs with rest length:
%tmp    = Xl(idxUp,:) - Xl;
%LUp    = sqrt(dot(tmp,tmp,2));
%TUp    = [(tmp(:,1) ./ LUp) (tmp(:,2) ./ LUp)];

%tmp    = Xl - Xl(idxDn,:);
%LDn    = sqrt(dot(tmp,tmp,2));
%TDn    = [(tmp(:,1) ./ LDn) (tmp(:,2) ./ LDn)];

%F      = zeros(numL, 2);
%F(:,1) = (K / (dtheta*dtheta)) * ( (LUp - L0Up) .* TUp(:,1) - (LDn - L0Dn) .* TDn(:,1) );
%F(:,2) = (K / (dtheta*dtheta)) * ( (LUp - L0Up) .* TUp(:,2) - (LDn - L0Dn) .* TDn(:,2) );


% target points
%F = F + Kp * ( XT - Xl );
