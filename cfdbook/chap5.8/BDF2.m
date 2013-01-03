% Stability domain of BDF2 scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 5.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.3  in the book

figure(1), clf

		% Plot of stability domain
theta = 0:0.01:1; theta = theta*pi;
z = 0.5*exp(-2*i*theta) - 2*exp(-i*theta) +3/2; plot(z)
