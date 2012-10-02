%% Example 2.4: Linear Transport and Diffusion
%
% In this example we will use operator splitting to numerically solve the
% equation:
% 
% $$u_t + \left(1+\sin(x)\right)u_x = \mu u_{xx}, \qquad x\in[-pi,pi]$$
% 
% That is, we split the evolution into two operators, transport and
% diffusion, that are solved for separately
%
% $$ S(t): \quad v_t + a(x) v_x = 0, \qquad v(x,0)=v_0(x)$$
%
% $$ H(t): \quad w_t = \mu w_{xx}, \qquad w(y,0)=w_0(y)$$
%
% Then the approximate solution is constructed from the formula
%
% $$u(x,t)\approx [H(\Delta t)\circ S(\Delta t) ]^n u_0(x)$$
%
% The two one-dimensional solution operators are approximated by the upwind
% method and a spectral difference method, respectively.

%% Initial setup
N  = 255;
x  = -pi+ 2*pi/N*(0:N);
u0 = 1.0*(abs(x)<1.0)+cos(x);
a  = 1+sin(x); 
T  = pi/2; 
%% Effect of splitting step
% In the first example, look at how the accuracy depends on the size of the
% splitting step. To this end, we compute the solution using 2, 8, 32, and
% 128 splitting steps and plot the solution at time t=0 (blue line), t=pi/4
% (green line), and time t=pi/2 (red line).
for i=1:4
   nstep=2*4^(i-1);
   u=transheat(u0,a,1.0,x,T,nstep);
   subplot(2,2,i); dplot=nstep/2;
   plot(x,u(:,1:dplot:nstep+1)); axis tight, title([num2str(nstep) ' steps']);
end
%%
% As we see from the plots, the solution is not very sensitive to the size
% of the splitting step and is captured quite accurately using only a few
% steps.

%% Width of shock layer
% The dynamics of this problem arises as a balance between advection and
% diffusion. This balance is controlled by the parameter \mu. In the next
% example, we look at different values of this parameter
mu = 10; nstep=32;
for i=1:4
   mu = 0.1*mu;
   u=transheat(u0,a,mu,x,T,nstep);
   subplot(2,2,i); dplot=nstep/2;
   plot(x,u(:,1:dplot:nstep+1)); axis tight, title(['mu = ' num2str(mu)]);
end
%%
% The plot shows how the transport effects become more and more dominant as
% \mu decreases toward zero. Notice, in particular the effect of the
% spatially dependent velocity for the two smaller \mu values.