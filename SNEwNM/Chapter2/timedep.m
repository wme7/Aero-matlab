% TIMEDEP This code solves the nonlinear parabolic pde
%
% u_t = u_xx + exp(u); u(0,t) = u(1,t) = 0; u(x,0) = 0; 0 < t < 1.
%
% with the backward Euler discretization. Newton's method is used
% for the nonlinear solver. The Jacobian is tridiagonal, so we
% use the banded differencing function.
%
% The value of u at the current time and the time step are passed
% to the nonlinear residual as MATLAB global variables.
%
% This problem is 1-D, so we can store the time history of the
% integration and draw a surface plot.
%
global uold dt
dt=.1;
nx=63; nt=1+1/dt;
dx=1/(nx+1);
tval=0:dt:1;
xval=0:dx:1;
%
% Use tight tolerances, Newton's method, and a tridiagonal Jacobian.
%
tol=[1.d-6,1.d-6];
parms=[40, 1, 0, 1, 1, 1];
uhist=zeros(nx+2,nt);
uold=zeros(nx,1); 
for it=1:nt-1
    [unew,it_hist,ierr]=nsold(uold,'ftime',tol,parms);
    uhist(2:nx+1,it+1)=unew;
    uold=unew;
end
%
% Plot the results.
%
mesh(tval,xval,uhist)
