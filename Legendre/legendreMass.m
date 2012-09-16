function M = legendreMass(k)
% Compute the mass matrix using legendre polinomials of order 'k' as a base
% for our expansion.

M = diag(2./(2*(0:k)+1));