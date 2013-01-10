% BVP2DEMO 
% This script solves the system of two-point boundary value 
% problems in Chapter 2 with nsold.m.
%
global L
L=20;
n=800;
u=zeros(n,1);
nh=n/2;
r=0:nh-1; h=L/(nh-1); r=r'*h;
%
% This choice of initial iterate gives the "correct" result.
% Try different initial iterates and
% watch Newton find a different solution! 
%
v=exp(-r.*r*.1); vp=-.2*r.*v;
u(1:2:n-1)=v; u(2:2:n)=vp;
tol=[1.d-12,1.d-12];
%
% Use Newton's method. The upper and lower bandwidths are both 2.
%
parms=[40, 1, 0, 1, 2, 2];
[sol, it_hist, ierr] = nsold(u,'bvpsys',tol,parms);
v=sol(1:2:n-1); vp=sol(2:2:n);
it_hist
plot(r,v,'-',r,vp,'--');
xlabel('t'); 
legend('v','v\prime');


