function ft=DGimexL(u)
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
global A b c Pleg w nx p dx dt u_alt gamma
pp=p+1;
v=reshape(u,nx,pp);
FS=zeros(nx,pp);
FC=zeros(pp,1);
ut=v;
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

            for i=1:nx
                FC(:)=v(i,:);
                FC=(FC'*Pleg).^2';
                for j=0:p
                    FS(i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2;
                end
            end
            %                 FC(:)=F(K,1,:);
            %             FU(:)=F(NV-K+1,1,:);
            %             FR(:)=FS(K,1,:);
            %             F_tmp(K,1,:)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c'-FR)' .* b);            
           ut(1,:)=( (FS(1,1:pp)')' .* b);
            for i=2:nx
                %                 FU(:)=F(K,i-1,:);
                %                 FC(:)=F(K,i,:);
                %                 FR(:)=FS(K,i,:);
                %                 F_tmp(K,i,:)=( (V(K)*A'*FC -V(K)* sum(FC) +V(K)* sum(FU) * c'-FR)' .* b);

                ut(i,:)=( (FS(i,1:pp)')' .* b);
            end

ft_a= v - u_alt - gamma*dt*ut;
ft=reshape(ft_a,nx*pp,1);