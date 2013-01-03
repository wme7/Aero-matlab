% Stability domain of Runge-Kutta method of Le and Moin JCP92:369-379

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.10  in the book

figure(1), clf, hold on

		% Plot of stability domain
global thetapar  center
theta = 0:0.01:0.7; theta = pi*theta; z = theta + i*theta;
center = -1/2; start =  abs(center);
for j = 1:length(theta)
  thetapar = theta(j); r = fzero('abspoly', start);
  start = r; z(j) = center + r*exp(i*theta(j));
end
plot(z)

		% Plot of imaginary axis
z = [0 3*i]; plot(z)
