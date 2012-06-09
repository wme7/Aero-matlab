function [delta] = evolve_kappa(phi, dx, dy, b, dx2, dy2)
%
% Finds the amount of evolution under a curvature-based
% force and using central differencing.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

phi_xp1 = zeros(size(phi));
phi_xm1 = zeros(size(phi));
phi_yp1 = zeros(size(phi));
phi_ym1 = zeros(size(phi));

phi_xp1(:,2:end) = phi(:,1:end-1);
phi_xp1(:,1) = 2*phi_xp1(:,2)-phi_xp1(:,3);
phi_xm1(:,1:end-1) = phi(:,2:end);
phi_xm1(:,end) = 2*phi_xm1(:,end-1)-phi_xm1(:,end-2);

phi_yp1(2:end,:) = phi(1:end-1,:);
phi_yp1(1,:) = 2*phi_yp1(2,:)-phi_yp1(3,:);
phi_ym1(1:end-1,:) = phi(2:end,:);
phi_ym1(end,:) = 2*phi_ym1(end-1,:)-phi_ym1(end-2,:);

phi_x = (phi_xm1 - phi_xp1)/(2*dx);
phi_y = (phi_ym1 - phi_yp1)/(2*dy);
phi_xx = (phi_xm1 - 2*phi + phi_xp1)/dx2;
phi_yy = (phi_ym1 - 2*phi + phi_yp1)/dy2;

dummy1 = circshift(phi, [1,1]);
dummy1(1,2:end) = 2*dummy1(2,2:end)-dummy1(3,2:end);
dummy1(:,1) = 2*dummy1(:,2)-dummy1(:,3);

dummy2 = circshift(phi, [-1,-1]);
dummy2(end,1:end-1) = 2*dummy2(end-1,1:end-1)-dummy2(end-2,1:end-1);
dummy2(:,end) = 2*dummy2(:,end-1)-dummy2(:,end-2);

dummy3 = circshift(phi, [1,-1]);
dummy3(1,1:end-1) = 2*dummy3(2,1:end-1)-dummy3(3,1:end-1);
dummy3(:,end) = 2*dummy3(:,end-1)-dummy3(:,end-2);

dummy4 = circshift(phi, [-1,1]);
dummy4(end,2:end) = 2*dummy4(end-1,2:end)-dummy4(end-2,2:end);
dummy4(:,1) = 2*dummy4(:,2)-dummy4(:,3);


phi_xy = (dummy1 + dummy2 - dummy3 - dummy4)/(4*dx*dy);

abs_grad_phi_sq = (phi_x.*phi_x + phi_y.*phi_y);
kappa_abs_phi = (phi_xx.*phi_y.*phi_y - 2.*phi_y.*phi_x.*phi_xy ...
    + phi_yy.*phi_x.*phi_x)./(abs_grad_phi_sq + (abs_grad_phi_sq == 0));

delta = b.*kappa_abs_phi;
