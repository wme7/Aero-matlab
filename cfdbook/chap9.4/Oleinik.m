% Application of Oleinik's entropy condition

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.4 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 9.11 in the book

		% ...............Input.........................................
a = 0.5;	% Parameter in flux function on page 371
		% ............End of input.....................................

x = 0.0:0.05:1.0; den = x.*x + 0.5*(1-x).^2; 
f = (x.*x)./den;	% Buckley-Leverett flux function

x1 = 1/sqrt(3); 
f1 = (x1*x1)/(x1*x1 + 0.5*(1-x1)^2);	%Point where tangent meets curve

figure(1), clf, hold on, plot(x,f,'-')
xx = [0  1];  xx = x1*xx; yy = [0  1]; yy = f1*yy; plot(xx,yy,'-')
zz = [0  f1]; xx = [x1 x1];			   plot(xx,zz,'--')
 
