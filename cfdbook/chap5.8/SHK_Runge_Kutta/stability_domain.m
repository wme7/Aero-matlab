% Stability domain of SHK Runge-Kutta scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Stability domain of SHK Runge-Kutta scheme

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.9  in the book
% Function called: abspoly

figure(1), clf, hold on
global thetapar  center

		% Plot of stability domain
theta = 0:0.01:1; theta = pi*theta; z = theta + i*theta;
center = -2; start =  abs(center);
for j = 1:length(theta)
  thetapar = theta(j); r = fzero('abspoly', start);
  start = r; z(j) = center + r*exp(i*theta(j));
end
plot(z)

cc = z(length(theta));		% Stability interval

		% Plot of imaginary axis
z = [0 3*i]; plot(z)

		% Plot of half ellipse
a = cc; b = 2.55; theta = -pi/2 +0.5*theta; z = a*cos(theta) - i*b*sin(theta);
plot(z,'--')

		% Plot of rectangle
z = [ i*c  (i*c-1)]; plot(z,'-.')
z = [ -1  (i*c-1)];  plot(z,'-.')
