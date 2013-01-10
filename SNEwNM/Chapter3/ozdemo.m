% OZDEMO
% This program creates the Ornstein-Zernike example in Chapter 3.
% [H,C]=OZDEMO returns the solution on a grid with a mesh
% spacing of 1/256.
%
function [h,c]=ozdemo
global L U rho
n=257;
epsilon=.1; sigma=2.0; rho=.2; beta=10; L=9;
dx=L/(n-1); r=0:dx:L; r=r'; 
% 
% Compute the potential and store it in a global variable.
%
U=elj(r,sigma,epsilon,beta);
%
tol=[1.d-8,1.d-8];
x=zeros(2*n,1);
parms=[40,80,-.1];
[sol, it_hist, ierr] = nsoli(x,'oz',tol);
%
% Unpack h and c.
%
h=sol(1:n); c=sol(n+1:2*n);
%
% Plot the solution.
%
figure(1)
subplot(1,2,1)
plot(r,h,'-');
ylabel('h','Rotation',1);
xlabel('r');
subplot(1,2,2)
plot(r,c,'-');
ylabel('c','Rotation',1);
xlabel('r');
%
% Do a second solve with constant forcing terms.
% Are you getting the same results each time? 
%
fb=it_hist(1,1);
parms=[40,80,-.1]; x=zeros(2*n,1);
[sola, it_hist1, ierr] = nsoli(x,'oz',tol,parms);
norm(sola-sol)
%
% plot residual vs iteration counter
%
figure(2)
ni=length(it_hist(:,1));
n1=length(it_hist1(:,1));
semilogy(0:ni-1,it_hist(:,1)/fb,'-',...
0:n1-1,it_hist1(:,1)/fb,'--');
legend('default','.1');
xlabel('Nonlinear iterations');
ylabel('Relative residual');
%
% plot residual vs function counter
%
figure(3)
semilogy(it_hist(:,2),it_hist(:,1)/fb,'-',...
it_hist1(:,2),it_hist1(:,1)/fb,'--');
xlabel('Function evaluations');
ylabel('Relative residual');
legend('default','.1');

%
function u=elj(r,sigma,epsilon,beta)
n2=length(r);
ra=r(2:n2);
r12=(sigma./ra).^12; r6=(sigma./ra).^6;
ua=exp(-4*beta*epsilon*(r12-r6));
u=[0,ua']';


