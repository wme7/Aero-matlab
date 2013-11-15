%% Example A.4: Linear Advection with Different Limiters
%
% The linear advection problem with periodic boundary conditions
%
% $$u_t + u_x = 0, \qquad u(x,0)=u_0(x), \qquad u(0,t)=u(1,t)$$
%
% is well suited for studying error mechanisms in numerical schemes for
% hyperbolic conservation laws. Here, we will use initial data consisting
% of a combination of a smooth, squared cosine wave and double step
% function to study the compressive and dissipative nature of four
% different limiters for a second-order nonoscillatory central-difference
% scheme.

%%
% If we choose the CFL number exactly equal to the stability limit of 0.5,
% the central scheme will produce solutions with excellent accuracy
% regardless of our choice of limiter. In practical computations, however,
% one cannot expect to simulate a linear wave with a CFL number equal the
% stability limit of the scheme. Hence, we choose a somewhat lower number
% to exhibit the typical behavior of the various limiters.
T   = 20;
CFL = 0.475;
N   = 100; h=1/N; x=h*(1:N);
u0  = (abs(x-.25)<=0.15).*(cos(pi*(10/3)*(x-0.25))).^2+...
   1.0*(abs(x-0.75)<0.15);
xx  = linspace(0,1,1001);
uf  = (abs(xx-.25)<=0.15).*(cos(pi*(10/3)*(xx-0.25))).^2+...
   1.0*(abs(xx-0.75)<0.15);

%%
% We consider four limiters (MinMod, vanLeer, MacCormack, and Superbee) and
% study the solution after 20 'passes' over the unit interval.
limiter = {'minmod', 'vanleer', 'mc', 'superbee'};
name   = {'MinMod', 'vanLeer', 'MacCormack', 'Superbee'};
for i=1:4
   u=central('linearfunc',u0,T,h,1,'periodic',CFL,limiter{i});
   subplot(2,2,i);
   plot(xx,uf,'-',x,u,'o','MarkerSize',4); axis([0 1 -0.2 1.3]);
   title(name{i});
end
%%
% The MinMod limiter is clearly dissipative, clipping the top of the smooth
% wave and smoothing the discontinuities of the double step, and in this
% way behaves somewhat like classical first-order scheme. The superbee
% limiter, on the other hand, picks steeper slopes and can thus resolve the
% discontinuities using very few cells but has also a tendency of
% overcompressing smooth linear waves, as observed for the smooth cosine
% profile. The other two limiters are somewhere in between.
