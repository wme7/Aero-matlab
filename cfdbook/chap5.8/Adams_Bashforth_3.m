% Stability domain of third order Adams-Bashforth scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.2  in the book

figure(1), clf, hold on

		% Plot of stability domain
theta = 0:0.01:1; theta = theta*pi;
z = exp(3*i*theta) - exp(2*i*theta);
z = 12*z./(23*exp(2*i*theta) - 16*exp(i*theta) + 5); plot(z)

		% Plot of half ellipse
a =6/11;b = (72/11)*sqrt(2/235); theta = pi/2 + 0.5*theta;
z = a*cos(theta) + i*b*sin(theta); plot(z,'--');

z = [0  0.8*i]; plot(z)% 	Plot of imaginary axis


