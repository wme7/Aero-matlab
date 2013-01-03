% Diagram of van Albada limiter

% Theory in Section 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Fig. 10.27 in the book

figure(1), clf, hold on

x = 0.0:0.01:4.0; y = (x.*x + x)./(1 + x.*x);	% Graph of van Albada limiter
plot(x,y)

x = 0.0:0.01:1.1; y = 0.5*(1+sqrt(2))*x; plot(x,y)
x = 0.0:0.01:4.0; y = 0.5*(1+sqrt(2))*ones(size(x)); plot(x,y)
y = ones(size(x)); plot(x,y)
x = 0.0:0.01:1.3; plot(x,x)
y = 0.0:0.01:1; x = ones(size(y)); plot(x,y,'--')
