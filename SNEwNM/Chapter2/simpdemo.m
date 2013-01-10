% SIMPDEMO
% This program solves the simple two dimensional problem in Chapter 2,
% and makes figure 2.1.
%
tol=[1.d-6,1.d-6];
%
% Create the mesh for the contour plot of \| f \|.
%
vl=.1:.5:2; vr=2:4:40; v=[vl,vr];
v=.5:4:40;
v=[.25,.5:2:40];
xr=-5:.2:5; n=length(xr); z=zeros(n,n);
for i=1:n
for j=1:n
    w=[xr(i),xr(j)]';
    z(i,j)=norm(simple(w));
end
end
%
% Newton's method
%
params=[40, 1, 0,0];
%
% x0 is a good initial iterate.
%
x0=[2,.5]';
[sn, errsn, ierrn, x_hist]=nsold(x0, 'simple', tol, params);
%
% x1 is a poor initial iterate. The iteration from x1 will stagnate
%    at a point where $F'$ is singular.
%
x1=[3,5]';
[sn2, errsn2, ierrn2, x_hist2]=nsold(x1, 'simple', tol, params);
%
% Draw a contour plot of $\| f \|$.
%
figure(1)
contour(xr,xr,z,v)
hold
%
% Use the x_hist array to plot the iterations on the coutour plot.
%
plot(x_hist(1,:),x_hist(2,:),'-*', x_hist2(1,:),x_hist2(2,:),'-o');
legend('Convergence','Stagnation');
xlabel('x_1');
ylabel('x_2');
axis([0 5 -5 5])
hold
%
% Default 
%
params=[40, 1000, .5, 0];
[sc, errsd, ierrd]=nsold(x0, 'simple', tol, params);
inewt=length(errsn(:,1)); cnewt=0:inewt-1;
idefault=length(errsd(:,1)); cdefault=0:idefault-1;
figure(2)
semilogy(cnewt,errsn(:,1),'-',cdefault,errsd(:,1),'--')
legend('Newton','default')
