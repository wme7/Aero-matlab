% HEQMOVIE This program creates a movie of the Broyden solution of
% the H-equation
%
function heqmovie
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
%
t1=cputime;
parms=[40,80];
[sol, it_hist, ierr,heq_hist] = brsola(x,'heq',tol,parms);
it_hist
%
% Build the movie of the iterations in h. Put the 
% iteration counter on the screen.
%
[hr,hc]=size(heq_hist)
for im=1:hc
    axis([0 1 0 3]);
    h=heq_hist(1:n,im);
    dx=1/(n-1); r=0:dx:1; r=r';
    plot(r,h,'-'); axis([0 1 0 3]); text(.5,.1,num2str(im));
    xlabel('x'); ylabel('h');
    mov(im)=getframe;
end
%
% Showtime!
%
movie(mov,1,1);
%
%
