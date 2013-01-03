% Stability domain of Adams-Bashforth_Euler scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.5  in the book

figure(1), clf, hold on

		% Plot of stability domain
theta = 0:0.01:1; theta = theta*pi;
z = - 2*(1-cos(theta)).^2./(3-cos(theta)) + i*2*sin(theta)./(3-cos(theta));
plot(z)

		% Stability domain of second order Adams-Bashforth scheme
z = exp(2*i*theta) - exp(i*theta); z = z./(1.5*exp(i*theta) - 0.5); plot(z,'-.')

		% Plot of oval
b = 1; a = 1; z =  - a*(1-cos(theta)) + i*b*(sin(theta)).^0.5; plot(z,'--')

		% Plot of half ellipse
a = 2.0; b = sqrt(1/3); theta = -pi/2 +0.5*theta;
z = - a*cos(theta) - i*b*sin(theta); plot(z,'.')
