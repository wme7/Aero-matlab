% Solution of Riemann problem for Buckley-Leverett equation

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.4 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 9.12 in the book

figure(1), clf, hold on
	% Left constant state
x = - 0.5:0.05:0; phi = ones(size(x)); plot(x,phi,'-')

	% Right constant state
x1 = 0.5 + 0.5*sqrt(3); x = x1 :0.05 :2; phi = zeros(size(x)); plot(x,phi,'-')

	% Expansion fan
phi = 0:0.05:1; phi = phi*(1 - 1/sqrt(3));  phi = 1/sqrt(3) + phi;
x = (phi - phi.^2)./((phi.^2 + 0.5*(1-phi).^2).^2); plot(x,phi,'-')
x = [x1   x1]; phi = [0   1/sqrt(3) ]; plot(x,phi,'-')
