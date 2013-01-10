% ATANEG MATLAB source for the arctangent function example.
%
% This example uses an initial iterate x0=10, which is far from
% the solution. Without a line search, Newton's method would
% diverge. newtsol has a line search and, as you can see by
% running the example, converges nicely.
%
% I use handle graphics to get the axes and fonts the way I
% want them. Comment out the "set" commands and look at the
% difference.
%
%
x0=10; tol=1.d-12;
[x,hist] = newtsol(x0, 'atan', tol,tol);
figure(1);
p1=subplot(1,1,1);
set(p1,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12],'FontSize',14);
semilogy(hist(:,1),abs(hist(:,2)),hist(2:5,1),abs(hist(2:5,2)),'o')
xlabel('Nonlinear iterations'); ylabel('Absolute Nonlinear Residual');
figure(2);
p2=subplot(1,1,1);
set(p2,'FontSize',14);
xran=-10:.1:10;
yran=atan(xran);
nz=length(xran);
xax=zeros(nz,1);
xvals=hist(:,4);
yvals=atan(xvals);
plot(xran,xax,'-',xran,yran,'-',xvals,yvals,'*');
%
% Here's a closeup view of convergence from a good initial iterate.
%
figure(3)
p3=subplot(1,1,1);
set(p3,'FontSize',14);
zran=-1.2:.01:1.2;
zvals=atan(zran);
nz=length(zran);
xaz=zeros(nz,1);
x0=1;
[x,hist] = newtsol(x0, 'atan', tol,tol);
nl=length(hist(:,4));
nx=2*nl-1;
xvals=zeros(1,nx);
yvals=zeros(1,nx);
xvals(1:2:nx)=hist(:,4);
xvals(2:2:nx)=hist(2:nl,4);
yvals(1:2:nx)=atan(hist(:,4));
plot(zran,xaz,'-',zran,zvals,'-',xvals,yvals,'-*');
xlabel('x');
ylabel('atan(x)');



