% Dissipation and dispersion plot for Lax-wendroff scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.3 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 9.1 in the book

		% ...............Input.........................................
sigma = [0.3  0.6  0.9]		% CFL numbers	
		% ..............End of input...................................
		
x = 0.0:0.05:1.0; theta = pi*x; dispersion = zeros(length(theta));

figure(1), clf
subplot(1,2,1), hold on, title('Dispersion','fontsize',18)
subplot(1,2,2), hold on, title('Dissipation','fontsize',18)
 
s = (sin(0.5*theta)).^2;
for k = 1:length(sigma)
  z = sigma(k)*sin(theta); 	z = z./(1-sigma(k)*sigma(k)*(1-cos(theta)));
  for j = 1:length(theta)
    if z(j) >= 0,    dispersion(j) = atan(z(j)) - sigma(k)*theta(j);
    else,	     dispersion(j) = atan(z(j)) + pi - sigma(k)*theta(j);
    end
  end
  dissipation = 1 - sqrt(1 + 4*sigma(k)^2*(sigma(k)^2 - 1)*s.^2);
  if k==1,lijntype='-';elseif k==2,lijntype = '--';else,lijntype='-.';end
  subplot(1,2,1), plot(x,dispersion,lijntype)
  subplot(1,2,2), plot(x,dissipation,lijntype)
end
