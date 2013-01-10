function fm=oz(x)
%OZ FM=OZ(X) is the nonlinear residual for the 
%   Ornstein-Zernike equations.
%
% Evaluate the nonlinearity as a sum of the substitution
% and compact parts.
%
fm=ozsub(x)+ozinteg(x);
%
% This is the substitution part. 
%
function km=ozsub(x)
global L U rho
n=length(x); n2=n/2; h=x(1:n2); c=x(n2+1:n);
nl=U.*exp(h-c);
mup=h+1- nl;
mdown=c+1-nl;
km=[mup',mdown']';
%
% This is the compact part of the nonlinearity. The upper
% component is zero, the lower is rho*L*c*h.
%
% rho and L are scalars for this problem
%
function km=ozinteg(x)
global L U rho
L3=L*L*L*rho;
n=length(x); n2=n/2; h=x(1:n2); c=x(n2+1:n);
kd=L3*convk(c,h);
kd(1)=2*kd(2)-kd(3);
kd(n2)=2*kd(n2-1)-kd(n2-2);
kup=zeros(n2,1);
km=[kup',kd']';
%
% Discrete convolution with the Hankel transform.
%
function hc=convk(h,c)
Lx=1.0;
n=length(h);
dr=Lx/n;
th=hank2(h); tc=hank2(c); fhc=th.*tc; hc=ihank2(fhc);
%
% Hankel transform using the fast sine transform.
%
function hf=hank2(f)
nf=length(f); n=nf-1; m=nf-2;
h=1/n; beta=2*(m+1)*(h^3); ff=beta;
p=1:m; p=p';
hft=ff*lsint(p.*f(2:n))./p;
hf=[0,hft',0]';
%
% Inverse Hankel transform using the fast sine transform.
%
function ihf=ihank2(f)
nf=length(f); n=nf-1; m=nf-2;
h=1/n; beta=2*(m+1)*(h^3); ff=2/(n*beta);
p=1:m; p=p'; 
hft=ff*lsint(p.*f(2:n))./p;
p1=1:n; p1=p1'; p1=p1.*p1; p1(n)=p1(n)/2;
ihf=[0,hft',0]';
% LSINT
% Fast sine transform with MATLAB's FFT.
%
function lf=lsint(f)
n=length(f);
ft=-fft([0,f']',2*n+2);
lf=imag(ft(2:n+1));

