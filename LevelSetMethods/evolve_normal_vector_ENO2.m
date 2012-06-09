function [delta, H1_abs, H2_abs] = evolve_normal_vector_ENO2(phi, dx, dy, Vn_ext, u_ext, v_ext)
%
% Finds the amount of evolution under a force in
% normal direction and a force based on a vector field,
% and using 2nd order accurate ENO scheme.
% Does not assume that phi is approximately a signed
% distance function and uses SLLF scheme.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


delta = zeros(size(phi)+4);
data_ext = zeros(size(phi)+4);
data_ext(3:end-2,3:end-2) = phi;

% Calculate the derivatives (both + and -)
phi_x_minus = zeros(size(phi)+4);
phi_x_plus = zeros(size(phi)+4);
phi_y_minus = zeros(size(phi)+4);
phi_y_plus = zeros(size(phi)+4);
% first scan the rows
for i=1:size(phi,1)
	phi_x_minus(i+2,:) = der_ENO2_minus(data_ext(i+2,:), dx);	
	phi_x_plus(i+2,:) = der_ENO2_plus(data_ext(i+2,:), dx);	
end

% then scan the columns
for j=1:size(phi,2)
	phi_y_minus(:,j+2) = der_ENO2_minus(data_ext(:,j+2), dy);	
	phi_y_plus(:,j+2) = der_ENO2_plus(data_ext(:,j+2), dy);	
end

[delta, H1_abs, H2_abs] = LLF_normal_vector(dx, dy, Vn_ext, u_ext, v_ext, phi_x_minus, phi_x_plus, phi_y_minus, phi_y_plus);
H1_abs = H1_abs(3:end-2,3:end-2);
H2_abs = H2_abs(3:end-2,3:end-2);
delta = delta(3:end-2,3:end-2);


