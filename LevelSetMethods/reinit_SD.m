function [phi] = reinit_SD(phi, dx, dy, alpha, accuracy, iterations)
%
% function [phi] = reinit_SD(phi, dx, dy, alpha, accuracy, iterations)
%
% Reinitializes phi into a signed distance function while preserving
% the zero level set (the interface or the curve).
%
% dx and dy are the resolution of the grid at x and y dimensions.
% alpha is a constant for calculating the euler step (dt). Should
% be between 0 and 1. 0.5 is quite safe whereas 0.9 can be risky.
% iterations specifies the number of iterations before the function returns.
% accuracy is the order of accuracy of derivative calculation.
% Allowed values for accuracy are 'ENO1', 'ENO2', 'ENO3', 'WENO'. 
% These correspond to 1st, 2nd, 3rd and 5th order accurate schemes 
% for calculating the derivative of phi.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


switch(accuracy)
	case 'ENO1'
		init_normal = @init_normal_ENO1;
		evolve_normal = @evolve_normal_ENO1;
	case 'ENO2'
		init_normal = @init_normal_ENO2;
		evolve_normal = @evolve_normal_ENO2;
	case 'ENO3'
		init_normal = @init_normal_ENO3;
		evolve_normal = @evolve_normal_ENO3;
	case 'WENO'
		init_normal = @init_normal_WENO;
		evolve_normal = @evolve_normal_WENO;
	otherwise
		error('Desired type of the accuracy is not correctly specified!');
end


S_phi_0 = phi./sqrt(phi.^2 + dx.^2);

Vn_ext = feval(init_normal, S_phi_0);
it=0;
t=0;
while(it < iterations)
	[delta_normal, H1_abs, H2_abs] = feval(evolve_normal, phi, dx, dy, Vn_ext);
	dt = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs);
	phi = phi + dt*(S_phi_0 - delta_normal);
	it = it+1;
	t = t+dt;
end


