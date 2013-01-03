function zL = extrap(z, limtype);	% Slope-limited extrapolation (MUSCL)
					% Staggered scheme
% Function called: limiter

% Numbering for pressure nodes:
%   1 1 2 2                     J J
% |-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|
% Numbering for momentum nodes:
% 1 1 2 2                        J-1  J
% |-o-|-o-|-o-|-o-|-o-|-o-|-o-|---o---|

J = length(z);
zL(1) = z(1); 
for j = 2:J-1
  b     = (z(j+1) - z(j)); a = (z(j) - z(j-1));
  zL(j) = z(j) + 0.5*limiter(a, b, limtype)*(z(j+1)-z(j));
end
zL(J) = z(J); 
