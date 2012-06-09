function [delta, H1_abs, H2_abs] = evolve_normal_vector_ENO1_SD(phi, dx, dy, Vn_ext, u_ext, v_ext)
%
% Finds the amount of evolution under a force in
% normal direction and a force based on a vector field,
% and using 1st order accurate ENO scheme.
% Assumes that phi is approximately a signed
% distance function and uses Godunov scheme.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

delta = zeros(size(phi)+2);
data_ext = zeros(size(phi)+2);
data_ext(2:end-1,2:end-1) = phi;

% Calculate the derivatives (both + and -)
phi_x_minus = zeros(size(phi)+2);
phi_x_plus = zeros(size(phi)+2);
phi_y_minus = zeros(size(phi)+2);
phi_y_plus = zeros(size(phi)+2);
phi_x = zeros(size(phi)+2);
phi_y = zeros(size(phi)+2);
% first scan the rows
for i=1:size(phi,1)
	phi_x_minus(i+1,:) = der_ENO1_minus(data_ext(i+1,:), dx);	
	phi_x_plus(i+1,:) = der_ENO1_plus(data_ext(i+1,:), dx);	
	phi_x(i+1,:) = select_der_normal_vector_SD(u_ext(i+1,:), Vn_ext(i+1,:), phi_x_minus(i+1,:), phi_x_plus(i+1,:));
end

% then scan the columns
for j=1:size(phi,2)
	phi_y_minus(:,j+1) = der_ENO1_minus(data_ext(:,j+1), dy);	
	phi_y_plus(:,j+1) = der_ENO1_plus(data_ext(:,j+1), dy);	
	phi_y(:,j+1) = select_der_normal_vector_SD(v_ext(:,j+1), Vn_ext(:,j+1), phi_y_minus(:,j+1), phi_y_plus(:,j+1));
end


abs_grad_phi = sqrt(phi_x.^2 + phi_y.^2);
H1_abs = abs(u_ext+Vn_ext.*phi_x./(abs_grad_phi + dx*dx*(abs_grad_phi<=1e-6)));
H2_abs = abs(v_ext+Vn_ext.*phi_x./(abs_grad_phi + dy*dy*(abs_grad_phi<=1e-6)));
H1_abs = H1_abs(2:end-1,2:end-1);
H2_abs = H2_abs(2:end-1,2:end-1);

delta = u_ext.*phi_x + v_ext.*phi_y + Vn_ext.*abs_grad_phi;
delta = delta(2:end-1,2:end-1);


