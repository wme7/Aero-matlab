% Dissipation and dispersion plot for kappa scheme

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.3 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Figs. 9.2 -- 9.4 in the book

		% ...............Input.........................................
kappa = [-1  0  1/3];			
sigma = 0.9;	% CFL number
		% ..............End of input...................................
		
x = 0.0:0.02:1.0; theta = pi*x; dispersion = x; s = (sin(theta./2)).^2;

figure(1), clf
subplot(1,2,1), hold on, title('Dispersion', 'fontsize',18)
subplot(1,2,2), hold on, title('Dissipation','fontsize',18)

for k = 1:length(kappa)
  a = sigma*(sigma - 1)*(1-kappa(k)+2*kappa(k)*sigma);
  znum=(sigma-a*s).*sin(theta); zden=(1-2*s*sigma^2+2*a*s.^2); z = znum./zden;
  for j = 1:length(theta)
    if z(j) >= 0,      dispersion(j) = atan(z(j));
    else,              dispersion(j) = (atan(z(j)) + pi);
    end
  end
  dispersion = dispersion - sigma*theta;  
  dissipation = sqrt(znum.^2 + zden.^2);  dissipation = 1 - dissipation;
  if k==1, ss ='--'; elseif k==2, ss ='-'; else, ss ='-.'; end  
  subplot(1,2,1),  plot(x,dispersion,ss)
  subplot(1,2,2),  plot(x,dissipation,ss)
end

