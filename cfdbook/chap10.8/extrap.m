function[zL,zR] = extrap(z, limtype);	% Slope-limited extrapolation
					% Function called: limiter 
J = length(z);
zL(1) = z(1); 
for j = 2:J-1
  b = (z(j+1) - z(j)); 	a = (z(j) - z(j-1));
  zL(j)   = z(j) + 0.5*limiter(a, b, limtype)*(z(j+1)-z(j));
  zR(j-1) = z(j) + 0.5*limiter(b, a, limtype)*(z(j-1)-z(j));
end
zR(J-1) = z(J);
