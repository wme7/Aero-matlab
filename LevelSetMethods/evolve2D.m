function [phi] = evolve2D(phi, dx, dy, alpha, iterations, accuracy, is_signed_distance, normal_evolve, Vn, vector_evolve, u, v, kappa_evolve, b)
%
% function [phi] = evolve2D(phi, dx, dy, iterations, accuracy, ...
%   is_signed_distance, normal_evolve, Vn, vector_evolve, u, v, kappa_evolve, b)
%
% Calculates evolution for a 2D curve (3D level set function) phi.
% phi is the input level set function. 
%
% dx and dy are the resolution of the grid at x and y dimensions.
% alpha is a constant for calculating the euler step (dt). Should
% be between 0 and 1. 0.5 is quite safe whereas 0.9 can be risky.
% iterations specifies the number of iterations before the function returns.
% normal_evolve, vector_evolve, kappa_evolve should be either set to 0 and 1. This
% indicates if these forces are present or not. If any of these are set to 1,
% corresponding  variables (e.g. Vn, u, v, b) should not be an empty array []. if
% either normal_evolve or vector_evolve are set to 1, accuracy needs to be specified,
% otherwise accuracy can be set to an empty array. Allowed values for accuracy are
% 'ENO1', 'ENO2', 'ENO3', 'WENO'. These correspond to 1st, 2nd, 3rd and 5th order
% accurate schemes for calculating the derivative of phi. if both normal_evolve and
% vector_evolve are set to 1, then is_signed_distance needs to be specified. If phi is
% approximately a signed distance function, set this variable to 1. If
% is_signed_distance == 1, Godunov scheme will be used, otherwise Stencil Local
% Lax-Friedrichs (SLLF) scheme will be used (See Osher&Fedkiw section 6.4).
%
% Other variables (these are either a scalar or of the same size as phi):
% Vn: Force in the normal direction.
% u: x component of the vector field
% v: y component of the vector field
% b: this specifies the weighting for the curvature (always positive everywhere).
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


tic

if alpha <= 0 | alpha >= 1 
    error('alpha needs to be between 0 and 1!');
end

% expand parameters if they are scalar
if length(Vn) == 1
	Vn = Vn.*ones(size(phi));
end
if length(u) == 1
	u = u.*ones(size(phi));
end
if length(v) == 1
	v = v.*ones(size(phi));
end
if length(b) == 1
	b = b.*ones(size(phi));
end

% Understand what kind of evolution the user is interested in.
if (normal_evolve == 1) & (vector_evolve == 1) & (kappa_evolve == 1)
	evolution_type = 'normal_vector_kappa';
	if isempty(Vn) | isempty(u) | isempty(v) | isempty(b)
		error('Some of the variables should not be empty!')
	end
	switch(accuracy)
		case 'ENO1'
			init_normal = @init_normal_ENO1;
			init_vector = @init_vector_ENO1;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO1_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO1;
			end
		case 'ENO2'
			init_normal = @init_normal_ENO2;
			init_vector = @init_vector_ENO2;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO2_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO2;
			end
		case 'ENO3'
			init_normal = @init_normal_ENO3;
			init_vector = @init_vector_ENO3;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO3_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO3;
			end
		case 'WENO'
			init_normal = @init_normal_WENO;
			init_vector = @init_vector_WENO;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_WENO_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_WENO;
			end
		otherwise
			error('Desired type of the accuracy is not correctly specified!');
	end
elseif (normal_evolve == 1) & (vector_evolve == 1) & (kappa_evolve == 0)
	evolution_type = 'normal_vector';
	if isempty(Vn) | isempty(u) | isempty(v)
		error('Some of the variables should not be empty!')
	end
	switch(accuracy)
		case 'ENO1'
			init_normal = @init_normal_ENO1;
			init_vector = @init_vector_ENO1;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO1_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO1;
			end
		case 'ENO2'
			init_normal = @init_normal_ENO2;
			init_vector = @init_vector_ENO2;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO2_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO2;
			end
		case 'ENO3'
			init_normal = @init_normal_ENO3;
			init_vector = @init_vector_ENO3;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_ENO3_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_ENO3;
			end
		case 'WENO'
			init_normal = @init_normal_WENO;
			init_vector = @init_vector_WENO;
			if is_signed_distance == 1
				evolve_normal_vector = @evolve_normal_vector_WENO_SD;
			else
				evolve_normal_vector = @evolve_normal_vector_WENO;
			end
		otherwise
			error('Desired type of the accuracy is not correctly specified!');
	end
elseif (normal_evolve == 1) & (vector_evolve == 0) & (kappa_evolve == 1)
	evolution_type = 'normal_kappa';
	if isempty(Vn) | isempty(b)
		error('Some of the variables should not be empty!')
	end
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
elseif (normal_evolve == 1) & (vector_evolve == 0) & (kappa_evolve == 0)
	evolution_type = 'normal';
	if isempty(Vn)
		error('Some of the variables should not be empty!')
	end
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
elseif (normal_evolve == 0) & (vector_evolve == 1) & (kappa_evolve == 1)
	evolution_type = 'vector_kappa';
	if isempty(u) | isempty(v) | isempty(b)
		error('Some of the variables should not be empty!')
	end
	switch(accuracy)
		case 'ENO1'
			init_vector = @init_vector_ENO1;
			evolve_vector = @evolve_vector_ENO1;
		case 'ENO2'
			init_vector = @init_vector_ENO2;
			evolve_vector = @evolve_vector_ENO2;
		case 'ENO3'
			init_vector = @init_vector_ENO3;
			evolve_vector = @evolve_vector_ENO3;
		case 'WENO'
			init_vector = @init_vector_WENO;
			evolve_vector = @evolve_vector_WENO;
		otherwise
			error('Desired type of the accuracy is not correctly specified!');
	end
elseif (normal_evolve == 0) & (vector_evolve == 1) & (kappa_evolve == 0)
	evolution_type = 'vector';
	if isempty(u) | isempty(v)
		error('Some of the variables should not be empty!')
	end
	switch(accuracy)
		case 'ENO1'
			init_vector = @init_vector_ENO1;
			evolve_vector = @evolve_vector_ENO1;
		case 'ENO2'
			init_vector = @init_vector_ENO2;
			evolve_vector = @evolve_vector_ENO2;
		case 'ENO3'
			init_vector = @init_vector_ENO3;
			evolve_vector = @evolve_vector_ENO3;
		case 'WENO'
			init_vector = @init_vector_WENO;
			evolve_vector = @evolve_vector_WENO;
		otherwise
			error('Desired type of the accuracy is not correctly specified!');
	end
elseif (normal_evolve == 0) & (vector_evolve == 0) & (kappa_evolve == 1)
	evolution_type = 'kappa';
	if isempty(b)
		error('Some of the variables should not be empty!')
	end
else
	error('Incorrect input parameters! normal_evolve, vector_evolve, kappa_evolve can only be 0 or 1 and all of them cannot be 0 at the same time');
end




% Now the main part of the program: Evolution
disp('Evolution type: ');
disp(evolution_type);
switch(evolution_type)
	case 'normal_vector_kappa'
		Vn_ext = feval(init_normal, Vn);
		[u_ext, v_ext] = feval(init_vector, u, v);
		[dx2, dy2] = init_kappa(dx, dy);
		it=0;
		t = 0;
		while(it < iterations)
			[delta_normal_vector, H1_abs, H2_abs] = feval(evolve_normal_vector, ...
				phi, dx, dy, Vn_ext, u_ext, v_ext);
			delta_kappa = evolve_kappa(phi, dx, dy, b, dx2, dy2);
			dt = get_dt_normal_vector_kappa(alpha, dx, dy, H1_abs, H2_abs, b, dx2, dy2);
			phi = phi + dt*(delta_kappa - delta_normal_vector);
			it = it+1
			t = t+dt;
		end
	case 'normal_vector'
		Vn_ext = feval(init_normal, Vn);
		[u_ext, v_ext] = feval(init_vector, u, v);
		it=0;
		t = 0;
		while(it < iterations)
			[delta_normal_vector, H1_abs, H2_abs] = feval(evolve_normal_vector, ...
				phi, dx, dy, Vn_ext, u_ext, v_ext);
			dt = get_dt_normal_vector(alpha, dx, dy, H1_abs, H2_abs);
			phi = phi - dt*delta_normal_vector;
			it = it+1
			t = t+dt;
		end
	case 'normal_kappa'
		Vn_ext = feval(init_normal, Vn);
		[dx2, dy2] = init_kappa(dx, dy);
		it=0;
		t = 0;
		while(it < iterations)
			[delta_normal, H1_abs, H2_abs] = feval(evolve_normal, phi, dx, dy, Vn_ext);
			delta_kappa = evolve_kappa(phi, dx, dy, b, dx2, dy2);
			dt = get_dt_normal_kappa(alpha, dx, dy, H1_abs, H2_abs, b, dx2, dy2);
			phi = phi + dt*(delta_kappa - delta_normal);
			it = it+1
			t = t+dt;
		end
	case 'normal'
		Vn_ext = feval(init_normal, Vn);
		it=0;
		t = 0;
		while(it < iterations)
			[delta_normal, H1_abs, H2_abs] = feval(evolve_normal, phi, dx, dy, Vn_ext);
			dt = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs);
			phi = phi - dt*delta_normal;
			it = it+1
			t = t+dt;
		end
	case 'vector_kappa'
		[u_ext, v_ext] = feval(init_vector, u, v);
		[dx2, dy2] = init_kappa(dx, dy);
		dt = get_dt_vector_kappa(alpha, dx, dy, u, v, b, dx2, dy2);
		it=0;
		t = 0;
		while(it < iterations)
			delta_vector = feval(evolve_vector, phi, dx, dy, u_ext, v_ext);
			delta_kappa = evolve_kappa(phi, dx, dy, b, dx2, dy2);
			phi = phi + dt*(delta_kappa - delta_vector);
			it = it+1
			t = t+dt;
		end
	case 'vector'
		[u_ext, v_ext] = feval(init_vector, u, v);
		dt = get_dt_vector(alpha, dx, dy, u, v);
		it=0;
		t = 0;
		while(it < iterations)
			delta_vector = feval(evolve_vector, phi, dx, dy, u_ext, v_ext);
			phi = phi - dt*delta_vector;
			it = it+1
			t = t+dt;
		end
	case 'kappa'
		[dx2, dy2] = init_kappa(dx, dy);
		dt = get_dt_kappa(alpha, dx, dy, b, dx2, dy2);
		it=0;
		t = 0;
		while(it < iterations)
			delta_kappa = evolve_kappa(phi, dx, dy, b, dx2, dy2);
			phi = phi + dt*delta_kappa;
			it = it+1;
			t = t+dt;
		end
	otherwise
		error('Something went wrong when checking the evolution type!');
end

t



toc












