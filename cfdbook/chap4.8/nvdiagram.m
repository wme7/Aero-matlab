% Normalized variable diagram

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 4.8 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 4.9  in the book

	% ......................Input..................................... 
kappa = [-1  0  1/3  1/2  1  1/3];
	% .....................End of input...............................

figure(1), clf
for j = 1:6		% Diagrams will be plotted for 6 schemes
  x = 0.0:0.01:1.0;
  y = (1 - kappa(j)/2)*x + (1+kappa(j))/4; % Characteristic of kappa-scheme
  subplot(2,3,j), hold on
  plot(x,y,'--'), axis([0.0  1.2  0.0 1.2])
  
  xx = 0.0:0.01:0.5; y = 2*xx;		   % Plot of dotted triangle
  plot(xx,y,':'), plot(x,x,':')
  xx = 0.5:0.01:1.0; y = ones(length(xx)); plot(xx,y,':'), plot(0.5,0.75,'o')
end

for j = 1:5  % Piecewise cubic-parabolic characteristic for general kappa scheme
  subplot(2,3,j), hold on, title(['kappa = ',num2str(kappa(j))])
  xx = 0.0:0.01:0.5; y = 2*xx +(kappa(j)-1)*xx.^2 - 2*kappa(j)*xx.^3;
  plot(xx,y)
  xx = 0.5:0.01:1.0; z = xx - ones(1,length(xx));
  if kappa(j) >= 0
    y = ones(1,length(xx)) + 0.5*kappa(j)*z + (kappa(j)-1)*z.^2;
  else
    y = ones(1,length(xx)) - (kappa(j)+1)*z.^2 - 2*kappa(j)*z.^3;
  end
  plot(xx,y)
end

subplot(2,3,6)		% Piecewise linear scheme
hold on, title(['PL'])
xx = 0.0:0.01:1.0; y=xx; kapfix = 1/3;
for j = 1:length(xx)
  f1 = max(min(2*xx(j),1),min(max(2*xx(j),1),xx(j)));
  f2 = (1-0.5*kapfix)*xx(j) + 0.25*(1+kapfix);
  y(j) = max(min(f1,f2),min(max(f1,f2),xx(j)));
end
plot(xx,y)


