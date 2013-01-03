function ye = profile(x)
a = - 0.3; b = -0.1; amplitude = 1;
if  x<= a |x>=b, y = 0;
else
 z = (x-a)/(b-a); y = amplitude*(1-(cos(pi*z))^2); 	% Smooth pulse
% y = 1;						% Block pulse
%  if z<= 0.5, y = 2*z; else, y = 1 - 2*(z-0.5); end	% Hat pulse
end
ye = y;
