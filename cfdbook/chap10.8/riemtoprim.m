function[u1,u2,u3] = riemtoprim(z1,z2,z3)
	% Transformation from Riemann invariants to primitive variables
global gamma
u = 0.5*(z1 + z3); 	c2 = (0.25*(gamma-1)*(z3 - z1)).^2;
u1 = (gamma*(exp(z2)./c2)).^(1/(1-gamma));	u2 = u.*u1;
u3 = u1.*(0.5*u.^2 + c2/(gamma*(gamma-1)));
