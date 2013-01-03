%NS1
% Numerical solution of 2-D instationary Navier-Stokes equations.
% Rectangular domain; nonuniform Cartesian staggered grid.
% Time discretization with omega-scheme.
% Pressure-correction method.
% Central or upwind scheme for inertia terms.

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Chapter 6 of
% 	Principles of Computational Fluid Dynamics, by P. Wesseling
% 	Springer-Verlag, Berlin etc., 2001. ISBN 3-540-67853-0
% See http://dutita0.twi.tudelft.nl/nw/users/wesseling/

% Theory is also given in Chapter 5 of
%	Computational Fluid Dynamics
%	Lecture Notes by P. Wesseling
%	Department of Applied Mathematical Analysis, ITS, TU Delft
% See http://ta.twi.tudelft.nl/nw/users/wesseling/

% This program generates Figures  and  in the Lecture Notes.

% Functions called: startup,problem_specification, grid_generation, draw_grid,
%		    viscous_matrix, grad_and_div, initial_condition, 
%		    right_hand_side, inertia_matrix, errornorm, isobarplot,
%		    streamlineplot, velocity_profile, relchangeplot  

startup, format compact, clear all
time0 = clock;
global geval n u0 v0 J K yu yv yseglen alpha
% ........................................Input..............................
geval = 3;	% Enter 1 for horizontal Poiseuille flow to the right
		%  or 2 for vertical Poiseuille flow
		%  or 3 for backward facing step
		%  or 4 for driven cavity
		%  or 5 for uniform flow under angle alpha without segmentation
		%  or 6 for uniform flow under angle alpha with segmentation
		%  or 7 for horizontal Poiseuille flow to the left
		
problem_specification	% Prepare file problem_specification.m as part of input
% User must prepare files ubd.m, vbd.m, pbd.m to specify boundary values 

% .................................End of input......................

grid_generation
draw_grid
h = waitbar(0,'Computing...');
viscous_matrix		% Generate viscous matrices Bu, Bv and Masku, Maskv
grad_and_div		% Generate pressure gradient matrices Pu, Pv and 
			%  divergence matrices Du, Dv
initial_condition	% Generate initial u0, v0, p0 
tic		 
[Lp,Up] = lu([Du Dv]*[Pu;Pv]);	% LU decomposition of pressure correction matrix
tijd = toc; disp(['Pres. corr. matrix LU time = ',num2str(tijd)])
if max(abs(sum(([Du Dv]*[Pu;Pv])'))) < 10*eps		
  % Do consistency check
end

nstep = floor(tend/dt);	% Number of time steps
%nstep = 2;
relchange = zeros(nstep,1);	% For plotting relative change per time step
t = 0; n = 0;
right_hand_side		% Generate right-hand sides ru, rv; has to be moved
			% inside time loop if boundary conditions time-dependent
u1 = u0; v1 = v0;
for n = 1:nstep		% Time-stepping with omega scheme
  t = n*dt;
  u0 = u1; v0 = v1;
%  right_hand_side	% Generate right-hand sides ru, rv; can be moved outside
			% time loop if boundary conditions not time-dependent

  % Navier-Stokes, predictor:
  inertia_matrix	% Generate inertia matrices and adapt right-hand sides 
  rum = (ru + rum)'; rvm = (rv + rvm)'; rum = rum(:); rvm = rvm(:);
  if n == 1, tic, end
  Cu = dt*(Bu + Cu); Cv = dt*(Bv + Cv);
  u1 = (speye(length(u0)) + omega*Cu)\...
        (Masku*u0 + dt*(rum - Pu*p0) - (1-omega)*Cu*u0);
  v1 = (speye(length(v0)) + omega*Cv)\...
        (Maskv*v0 + dt*(rvm - Pv*p0) - (1-omega)*Cv*v0);
  
  % End of Navier-Stokes predictor
   
  % Stokes:
  %rum = ru'; rum = rum(:); rvm = rv'; rvm = rvm(:);  
  %if n == 1, tic, end 
  %u1 = Bu\(dt*(rum - Pu*p0) + Masku*u0); 
  %v1 = Bv\(dt*(rvm - Pv*p0) + Maskv*v0);
  % End of Stokes
  
  if n == 1, tijd = toc; disp(['velocity solve time = ',num2str(tijd)]), tic,end
  dp = Up\(Lp\(Du*u1 + Dv*v1));	p0 = p0 + dp/dt;
  u1 = u1 - Pu*dp;	v1 = v1 - Pv*dp;
  if n == 1, tijd = toc; disp(['press. corr. time = ',num2str(tijd)]),end
  mm = max([norm(u1,inf),norm(v1,inf)]); relchange(n) = max([norm(u1-u0,inf)...
  /mm,norm(v1-v0,inf)/mm, norm(dp,inf)/(norm(p0,inf) + 0.5*mm^2)]);
  waitbar(n/nstep)
end

disp(['Relative change = ',num2str(relchange(end))])
divergence = norm(Du*u1 + Dv*v1,inf);
if (geval==1)|(geval==2)
  errornorm		% Maximum norm of error in outflow profile if exact
			% solution is available
end
close(h)
tic
isobarplot		% Plot of isobars and velocity vectors
streamlineplot 		% Plot of streamlines and velocity vectors
velocity_profile	% Plot of velocity profiles
relchangeplot		% Plot of relative change per time step
tijd = toc; disp(['Plotting time = ',num2str(tijd)])
disp(['Total time = ',num2str(etime(clock,time0))])
