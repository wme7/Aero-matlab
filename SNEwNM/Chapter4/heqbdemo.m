% HEQBDEMO Solve the H-equation with brsola
% This program creates the H-equation example in Chapter 4.
%
%
global A_heq;
c=.9;
n=100;
%
% nodal points for the midpoint rule
%
gr=1:n;
gr=(gr-.5)/n;
gr=gr';
%
% form and store the kernel of the integral operator
%
cc=.5*c/n;
A_heq=ones(n,1)*gr'; A_heq=cc*A_heq'./(A_heq+A_heq');
%
tol=[1.d-8,1.d-8];
x=ones(n,1);
[sol, it_hist, ierr] = brsola(x,'heq',tol);
nb=length(it_hist(:,1));
[sol, it_histn, ierr] = nsoli(x,'heq',tol);
nn=length(it_histn(:,1));
[sol, it_histd, ierr] = nsold(x,'heq',tol);
nd=length(it_histd(:,1));
semilogy(0:nb-1,it_hist(:,1)/it_hist(1,1),'-',...
0:nn-1,it_histn(:,1)/it_histn(1,1),'--',...
0:nd-1,it_histd(:,1)/it_histd(1,1),'-.');
legend('brsola','nsoli','nsold');
xlabel('Nonlinear iterations');
ylabel('Relative residual');

