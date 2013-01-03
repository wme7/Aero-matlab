% Stability domain of second order extrapolated BDF scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.7  in the book

figure(1), clf, hold on

		% Plot of stability domain
theta = 0:0.01:0.5; theta = theta*pi;
a = (6*cos(theta)-9)./(2*cos(2*theta) - 4*cos(theta)) - 3/2;
b = (sin(theta).*(2-cos(theta)))./(cos(2*theta) - 2*cos(theta));
z = - a - i*b; plot(z)

		% Plot of oval
a = 1; b =(2*a/4)^(1/4); theta = 0:0.002:1; theta = theta*pi;
z =  - a*(1-cos(theta)) + i*b*(sin(theta)).^0.5; plot(z,'--')

		% Plot of parabola
yy = 0:0.02:1.8; b = 1; par = -(yy/b).^2; plot(par,yy,'.')
