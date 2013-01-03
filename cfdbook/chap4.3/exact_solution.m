% Exact solution of 1-D convection-diffusion equation with Dirichlet 
% boundary conditions:
%	du/dx - (1/Pe)d2u/dx2 = 0,   u(0) = a,   u(1) = b 
function res = exact_solution(x,pe,a,b)
res = a + (b-a)*(exp((x-1)*pe) - exp(-pe))/(1 - exp(-pe));
