function ft=DGimex(u)
% PDETIME 
% FT = PDETIME(U) is the nonlinear residual for time-dependent 
% problem in Chapter 3.
% This code has the zero boundary conditions, the right
% side, and the preconditioner built in.
%
% The equation is u_t = (u_xx + u_yy) - 20*u*(u_x + u_y)  + f
%
% Where f is the right hand side of the steady-state example
%
%
global b Pleg w p dx dt u_alt gamma i
pp=p+1;
v=u;
FS=zeros(1,pp);
FC=zeros(pp,1);
% ut=v;
%
% Left preconditioned nonlinear residual for implicit Euler discretization.
%
% v=20*u.*(dxmf(u)+dymf(u))-rhsf;
% %
% %
% %ft=dt*u + fish2d(u - uold + dt*v);
% ft=dt*u + fish2d(u - u_alt + dt*v);

%
%
%ft=dt*u + fish2d(u - uold + dt*v);
%ft= u - u_alt + gamma*dt*(v+lapmf(u));

                FC(:)=v(:);
                FC=(FC'*Pleg).^2';
                for j=0:p
                    FS(1,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2;
                end
                ut=( (FS(1,1:pp)')' .* b)';
                ut= FS(1,1:pp)' .* b';
v=u_alt(i,:)';
ft= u - v - gamma*dt*ut;

