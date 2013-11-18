%% Project: MAC Scheme for Stokes Equations
%
% The purpose of this project is to implement the simple and popular MAC
% scheme for solving Stokes equations in two dimensions.
%
% We implement matrix free version of MAC scheme on rectangular grids. The
% unknowns are represented by matrices; see <http://math.uci.edu/~chenlong/226/MACStokes.pdf MAC for Stokes Equations>.

%% Step 1: Gauss-Seidel relaxation of velocity
% Given a pressure approximation, relax the momentum equation to update
% velocity.
% See
% <http://www.math.uci.edu/~chenlong/iFEM/doc/project/html/projectMG.html Multigrid> project on the matrix free implemenation of G-S relaxation.
% Note that the boundary or near boundary dof should be updated diffeently.
% The stencil should be changed according to different boundary conditions.

%% Step 2: Distributive relaxation of velocity and pressue 
%
% # form the residual for the continuity equation: rc = div u.
% # solve the Poisson equation for pressure Ap*dp = rc by one G-S.
% # distribute the correction to velocity by u = u + grad dp;
% # update the pressure by p = p - Ap*dp;
%
% Every step can be implemented in a matrix-free version; see
% <http://math.uci.edu/~chenlong/226/MGStokes.pdf Multigrid Methods for Stokes Equations>
%
% An alternative approach is: for each cell
%
% # form the residual for the continuity equation: rc = div u.
% # update velocity at four edges such that divu = 0 in this cell.
% # update pressure in the 5-cells to keep the residual of momentum
% equation unchanged.
%
% Exact formulae can be found in
% <http://math.uci.edu/~chenlong/226/MGStokes.pdf Multigrid Methods for Stokes Equations>

%% Step 3: Test Example
%
% We use a simple model of colliding flow with analytic solutions to test
% the code. The domain is [-1,1]^2 and the analytic solution is:
%
% $$u = 20xy^3; v = 5x^4 - 5y^4; p = 60x^2y - 20y^3 + constant.$$
%
% Compute the data f and Dirichlet boundary condition g_D and solve Stokes
% equation on the unit square using methods in Part I,II,II and check the
% rate of convergence.