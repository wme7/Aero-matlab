% HEQKDEMO Solve the H-equation with nsoli.
% This program creates the H-equation example in Chapter 3.
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
%
% GMRES
%
x=ones(n,1); parms = [40,40,.9,1];
[sol, it_histg, ierr] = nsoli(x,'heq',tol,parms);
%
% BICGSTAB
%
x=ones(n,1); parms = [40,40,.9,3];
[sol, it_histb, ierr] = nsoli(x,'heq',tol,parms);
%
% TFQMR
%
x=ones(n,1); parms = [40,40,.9,4];
[sol, it_histt, ierr] = nsoli(x,'heq',tol,parms);
figure(1)
ng=length(it_histg(:,1));
nb=length(it_histb(:,1));
nt=length(it_histt(:,1));
semilogy(0:ng-1,it_histg(:,1)/it_histg(1,1),'-',...
0:nb-1,it_histb(:,1)/it_histg(1,1),'--',...
0:nt-1,it_histt(:,1)/it_histg(1,1),'-.');
xlabel('Nonlinear iterations');
ylabel('Relative residual norm');
legend('GMRES','BICGSTAB','TFQMR');
figure(2)
semilogy(it_histg(:,2),it_histg(:,1)/it_histg(1,1),'-',...
it_histb(:,2),it_histb(:,1)/it_histb(1,1),'--',...
it_histt(:,2),it_histt(:,1)/it_histt(1,1),'-.');
xlabel('Function evaluations');
ylabel('Relative residual norm');
legend('GMRES','BICGSTAB','TFQMR');

