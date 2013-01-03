function[z1,z2,z3] = primtoriem(u1,u2,u3)
	% Transformation from primitive to Riemann invariants
global gamma
u = u2./u1; 			       c2 = gamma*(gamma-1)*(u3./u1 - 0.5*u.^2);
z2 = log((1/gamma)*c2.*u1.^(1-gamma)); c2 = sqrt(c2);
z1 = u - 2*c2/(gamma-1);	       z3 = u + 2*c2/(gamma-1);
