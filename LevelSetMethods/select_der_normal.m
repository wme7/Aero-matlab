function [der] = select_der_normal(Vn, der_minus, der_plus)
%
% Under a force in the normal direction,
% select a derivative value given (plus) and (minus) derivatives.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

if size(der_minus) ~= size(der_plus) | size(der_plus) ~= size(Vn)
    error('plus, minus derivative vectors and normal force (Vn) need to be of equal length!');
end

der = zeros(size(der_plus));

for i=1:numel(Vn)
	Vn_der_m = Vn(i)*der_minus(i);
	Vn_der_p = Vn(i)*der_plus(i);
	if Vn_der_m <= 0 & Vn_der_p <= 0
		der(i) = der_plus(i);
	elseif Vn_der_m >= 0 & Vn_der_p >= 0
		der(i) = der_minus(i);
	elseif Vn_der_m <= 0 & Vn_der_p >= 0
		der(i) = 0;
	elseif Vn_der_m >= 0 & Vn_der_p <= 0
		if abs(Vn_der_p) >= abs(Vn_der_m)
			der(i) = der_plus(i);
		else
			der(i) = der_minus(i);
		end
	end
end





