% Graph of Buckley-Leverett flux function

% Theory in Section 9.4 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 9.10 in the book

		% ...............Input.........................................
a = 0.5;	% Parameter in flux function on page 371
		% ............End of input.....................................
x = 0.0:0.05:1.0; den = x.*x + 0.5*(1-x).^2; f = (x.*x)./den;
figure(1), clf, plot(x,f,'-')
