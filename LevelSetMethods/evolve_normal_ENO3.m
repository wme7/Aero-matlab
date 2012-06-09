function [delta, H1_abs, H2_abs] = evolve_normal_ENO3(phi, dx, dy, Vn)
%
% Finds the amount of evolution under a force in
% normal direction and using 3rd order accurate ENO scheme
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

delta = zeros(size(phi)+6);
data_ext = zeros(size(phi)+6);
data_ext(4:end-3,4:end-3) = phi;

% Calculate the derivatives (both + and -)
phi_x_minus = zeros(size(phi)+6);
phi_x_plus = zeros(size(phi)+6);
phi_y_minus = zeros(size(phi)+6);
phi_y_plus = zeros(size(phi)+6);
phi_x = zeros(size(phi)+6);
phi_y = zeros(size(phi)+6);
% first scan the rows
for i=1:size(phi,1)
	phi_x_minus(i+3,:) = der_ENO3_minus(data_ext(i+3,:), dx);	
	phi_x_plus(i+3,:) = der_ENO3_plus(data_ext(i+3,:), dx);	
	phi_x(i+3,:) = select_der_normal(Vn(i+3,:), phi_x_minus(i+3,:), phi_x_plus(i+3,:));
end

% then scan the columns
for j=1:size(phi,2)
	phi_y_minus(:,j+3) = der_ENO3_minus(data_ext(:,j+3), dy);	
	phi_y_plus(:,j+3) = der_ENO3_plus(data_ext(:,j+3), dy);	
	phi_y(:,j+3) = select_der_normal(Vn(:,j+3), phi_y_minus(:,j+3), phi_y_plus(:,j+3));
end

abs_grad_phi = sqrt(phi_x.^2 + phi_y.^2);

H1_abs = abs(Vn.*phi_x.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H2_abs = abs(Vn.*phi_y.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H1_abs = H1_abs(4:end-3,4:end-3);
H2_abs = H2_abs(4:end-3,4:end-3);

delta = Vn.*abs_grad_phi;
delta = delta(4:end-3,4:end-3);


