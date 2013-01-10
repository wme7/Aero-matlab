% OZMOVIE  This program creates a movie of the Newton-Krylov solution of
% the Ornstein-Zernike equations.
%
function ozmovie
global L U rho
n=257;
epsilon=.1; sigma=2.0; rho=.2; beta=10; L=9;
dx=L/(n-1); r=0:dx:L; r=r'; 
%
% U is a global variable used in the function evaluation. See the text.
%
U=elj(r,sigma,epsilon,beta);
%
tol=[1.d-8,1.d-8];
x=zeros(2*n,1);
t1=cputime;
parms=[40,80,-.1];
[sol, it_hist, ierr,oz_hist] = nsoli(x,'oz',tol,parms);
%
% Build the movie of the iterations in h. Put the 
% iteration counter on the screen.
%
[hr,hc]=size(oz_hist); n=hr/2;
for im=1:hc
    axis([0 9 -4 4]);
    h=oz_hist(1:n,im);
    L=9; dx=L/(n-1); r=0:dx:L; r=r';
    plot(r,h,'-'); axis([0 9 -2 3]); text(3,-1,num2str(im));
    xlabel('r'); ylabel('h');
    mov(im)=getframe;
end
%
% Showtime!
%
movie(mov,1,1);
%
%
%
function u=elj(r,sigma,epsilon,beta)
n2=length(r);
ra=r(2:n2);
r12=(sigma./ra).^12; r6=(sigma./ra).^6;
ua=exp(-4*beta*epsilon*(r12-r6));
u=[0,ua']';
