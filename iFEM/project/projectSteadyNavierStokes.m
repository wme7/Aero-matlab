%% Project: Steady Navier-Stokes Equations
%
% The purpose of this project is to implement multigrid methods for solving
% steady state Navier-Stokes equations in two dimensions.

%% Problem Setting: Driven Cavity
%
% Driven cavity problem. The domain is [-1,1]^2. Navier-Stokes equation
% with zero source and zero Dirichlet boundary condition except on the top:
%
% $$ \{ y=1; 0<= x <= 1 | u = 1, v = 0 \} $$

%% Step 1: Multigrid of MAC for Stokes Equations
%
% This is the Part III of Project: <projectMGStokes.html Fast Solvers for
% Stokes Equations>. You can download my implementation:
% <http://math.uci.edu/~chenlong/code/StokesMAC2012spring.zip StokesMAC.zip>
% Run squareStokesMAC.m to check the almost second order convergence. Read
% StokesVcycle.m to understand components of the algorithm.

%% Step 2: FAS for Nonlinear Convection Diffusion Equation
%
% Implement FAS for solving the two-point boundary value problem
%
% $$ - u''(x) + u(x)u'(x) = f(x), 0<x<1, u(0) = u(1) = 0. $$
%
% The nonlinear Gauss-Seideal using the centeral difference scheme can be
% found in the book: A Multigrid Tutorial.

%% Step 3: FAS for Navier-Stokes Equations with low Reynold Number
%
% Combine code from Step 1 and Step 2 to solve the Driven Cavity problem
% with low Reynold number or equivalently big visicosity constant.
%
% Test your FAS for nu = 1 first and then try nu = 4h, 2h, h.

%% Step 4: Defect correction for High Reynold Number
%
% Try your FAS for nu = 1e-6. If you want to keep working on this, please
% talk to me.