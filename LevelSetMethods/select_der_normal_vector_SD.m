function [der] = select_der_normal_vector_SD(u, Vn, der_minus, der_plus)
%
% Under a force in the normal direction and a force based on a 
% vector field, and assuming phi is approximately a signed distance
% function, select a derivative value given (plus) and (minus) derivatives.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

if size(der_minus) ~= size(der_plus) | size(der_plus) ~= size(Vn) | size(Vn) ~= size(u)
    error('plus, minus derivative vectors, vector field component (u) and normal force (Vn) need to be of equal length!');
end

der = zeros(size(der_plus));

H1_m = u + Vn.*der_minus;
H1_p = u + Vn.*der_plus;
absH1_m = abs(H1_m);
absH1_p = abs(H1_p);
for i=1:numel(Vn)
	if H1_m(i) >= 0 & H1_p(i) >= 0
		der(i) = der_minus(i);
	elseif H1_m(i) <= 0 & H1_p(i) <= 0
		der(i) = der_plus(i);
	elseif H1_m(i) <= 0 & H1_p(i) >= 0
		der(i) = -u(i)/Vn(i);
	elseif H1_m(i) >= 0 & H1_p(i) <= 0
		if absH1_m(i) >= absH1_p(i)
			der(i) = der_minus(i);
		else
			der(i) = der_plus(i);
		end
	end
end





