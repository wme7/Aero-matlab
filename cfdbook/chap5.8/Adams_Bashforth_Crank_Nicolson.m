% Stability domain of Adams-Bashforth-Crank-Nicolson scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.6  in the book

figure(1), clf, hold on

		% Plot of stability domain
theta = 0:0.005:0.8; theta = theta*pi; c = (2-cos(theta)).*(1+cos(theta));
a = 0.5*(1-cos(theta)).^2./c; b = - 2*sin(theta)./c; z = -(a+i*b); plot(z)

		% Plot of oval
a = 0.5; b = 0.75^0.25; theta = 0:0.002:1; theta = theta*pi;
z =  - a*(1-cos(theta)) + i*b*(sin(theta)).^0.5; plot(z,'--')

		% Plot of parabola
yy = 0:0.02:2; par = -0.75*(yy).^2; plot(par,yy,'.')
