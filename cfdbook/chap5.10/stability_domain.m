% Stability domain of multistep omega scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 5.10 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Amplification factor is giveb by equation (5.103)

omega = 1 - sqrt(0.5);			% Parameters in multistage omega scheme
alpha = (1 - 2*omega)/(1 - omega);
alphap = 1 - alpha;
omegap = 1 - 2*omega;

x = -0.03:0.01:0.5; y = 0:0.5:9; [X,Y] = meshgrid(x,y); Z = X + i*Y;

	% Absolute value of g(z)
One = ones(size(Z));
F = abs(((One + alphap*omega*Z).^2).*(One + alpha*omegap*Z)...
       .*((One - alpha*omega*Z).^(-2)).*((One - alphap*omegap*Z).^(-1)));
       
figure(1), clf, hold on
v = [ 0.9 0.95 1.0 ];   c = contour(x,y,F,v);   clabel(c,'fontsize',18);
line([0 0],[0 9])		% Draw imaginary axis
title('Contour plot of |g(z|','fontsize',18)

	% Value of abs[g(z)] on imaginary axis
z = i*3;    ff = abs(((1 + alphap*omega*z)^2)*(1 + alpha*omegap*z)...
      *((1 - alpha*omega*z)^(-2))*((1 - alphap*omegap*z)^(-1)))
