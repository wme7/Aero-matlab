%% Example A.7: Scalar Front Tracking
% In this example we show the basic concept of front tracking. The main
% idea of this method is to take a conservation law
%
% $$ u_t + f(u)_x = 0, \qquad u(x,0)=u_0(x)$$
%
% and replace the flux function "f" by a piecewise linear approximation "g"
% and the initial data "u0" by a piecewise constant approximation "v0".
% Then the new conservation law
%
% $$ v_t + g(v)_x = 0, \qquad v(x,0)=v_0(x)$$
%
% can be solved analytically using an algorithm with a finite number of
% steps if "v0" has a finite number of discontinuities. When used as a
% numerical method, the front tracking is unconditionally stable in the
% sense that there is no restriction on the time step. Moreover, the method
% has first-order convergence also for discontinuous solutions.
%
% As an example problem, we consider a cubic flux function and a sine wave
% truncated to the interval [-1,1].

%% Definition of domain and time
T = 1;
umin=-1; umax=1.0;
xmin=-2; xmax=2;
addpath Scalar_Fronttracking;

%% Approximation of flux function
% The accuracy of the front tracking method is defined two parameters: the
% parameter determining how well we approximate the flux function and the
% parameter that determines how well we approximate the initial function.
% The flux approximation is determined by the number of breakpoints
% "nbreak", which is equal the number of piecewise linear segments plus
% one.
nbreak=15;
uf=linspace(umin,umax,nbreak);
v = linspace(umin,umax,513);
ff=uf.^3;
subplot(1,2,1)
plot(v, v.^3, '-', uf,ff,'-o','MarkerSize',4);
xlabel('u'); ylabel('f(u)');  title(' Flux function ');
subplot(1,2,2)
plot(v, v.^3, '-', uf,ff,'-o','MarkerSize',4), axis([-0.4 0.4 -0.05 0.05])
xlabel('u'); ylabel('f(u)');  title(' Zoom');


%% Approximation of initial data
% To approximate the initial data, we choose initial values among the
% breakpoints of the flux function; that is, we approximate u0 along the
% "u-axis" rather than along the "x-axis". Doing so, the solution of the
% Riemann problem can (in principle) be performed in interger arithmetic
% and thus be *very* fast.
xx = linspace(xmin,xmax,513);
clf, plot(xx,initial(xx),'-r')

ndisc=100;
[u0,x0]=initialvalues(uf,xmin,xmax,'initial',ndisc);
plotFTresult(u0,x0);

%% Putting it all together
% We plot the approximation of the flux function and the initial data. The
% we perform front tracking: (i) solve all Riemann problems defined by
% discontinuities in the initial data, (ii) track all evolving fronts
% defined by the piecewise linear flux segments that are part of each local
% Riemann solution and compute potential collision, (iii) whenever two such
% fronts collide, solve a new Riemann problem. Repeat steps (ii) and (iii)
% until there are no more collisions. Our implementation plots the fronts
% while tracking them. Finally, we plot the analytical, piecewise constant
% solution of the perturbed PDE problem

% Piecewise linear flux
subplot(2,2,1), plot(uf,ff,'-o','MarkerSize',4)
xlabel('u'), ylabel('f(u)'), title(' Flux function ')

% Piecewise constant initial function
subplot(2,2,2), plotFTresult(u0,x0)
xlabel('x'), ylabel('u_0(x)'), title(' Initial data ')

% Tracking the fronts and plotting while doing this.
subplot(2,2,3), [u,x]=SFrontTrack(u0,x0,T,uf,ff);
xlabel('x'), ylabel('t'), title(' Fronts in the \it (x,t) \rm plane')

% Plot the exact solution of the perturbed problem
subplot(2,2,4), plotFTresult(u,x) 
xlabel('x'); ylabel(strcat('u(x,',num2str(T),')')); 
title(' Solution ');

%% Increasing the resolution
% We now repeate the same excercise with increased resolution

% Piecewise linear flux
clf, subplot(2,2,1), 
uf=linspace(umin,umax,101); ff=uf.^3;
plot(uf,ff,'-o','MarkerSize',4)
xlabel('u'), ylabel('f(u)'), title(' Flux function ')

% Piecewise constant initial function
subplot(2,2,2), 
[u0,x0]=initialvalues(uf,xmin,xmax,'initial',200);
plotFTresult(u0,x0)
xlabel('x'), ylabel('u_0(x)'), title(' Initial data ')

% Tracking the fronts and plotting while doing this.
subplot(2,2,3), [u,x]=SFrontTrack(u0,x0,T,uf,ff);
xlabel('x'), ylabel('t'), title(' Fronts in the \it (x,t) \rm plane')

% Plot the exact solution of the perturbed problem
subplot(2,2,4), plotFTresult(u,x) 
xlabel('x'); ylabel(strcat('u(x,',num2str(T),')')); 
title(' Solution ');


%% Disclaimer
% Herein, we have made a simple implementation in MATLAB to demonstrate the
% principles behind front tracking. The implementation is made using
% manipulation of arrays, which is quite slow, and the computational speed
% of the implementation does not reflect the high efficiency one can
% observe for a front-tracking method implemented using pointers in C or
% C++.