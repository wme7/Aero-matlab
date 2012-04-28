%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   2D lid driven cavity by using LBE Multi-Relaxation time scheme    %%%
%%%   Writer : Arman Safdari                                            %%%
%%%   E-mail : sarman2@live.utm.my                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
%%%%  definitions %%%%%%%
lx=100;
ly=100;
rho(1:lx,1:ly)=1;
Uo=0.1;
Nu=0.01;
Re=abs(Uo)*lx/Nu;
tau=3*Nu+0.5;

u=zeros(lx,ly);
v=zeros(lx,ly);
mEq=ones(9,lx,ly);
f=ones(9,lx,ly);

ex=[0,1,0,-1,0,1,-1,-1,1];
ey=[0,0,1,0,-1,1,1,-1,-1];
S=[1,1.4,1.4,1,1.2,1,1.2,1/tau,1/tau];
M=[1, 1, 1, 1, 1, 1, 1, 1, 1;-4,-1,-1,-1,-1, 2, 2, 2, 2; 
   4,-2,-2,-2,-2, 1, 1, 1, 1; 0, 1, 0,-1, 0, 1,-1,-1, 1;
   0,-2, 0, 2, 0, 1,-1,-1, 1; 0, 0, 1, 0,-1, 1, 1,-1,-1;
   0, 0,-2, 0, 2, 1, 1,-1,-1; 0, 1,-1, 1,-1, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 1,-1, 1,-1];
InvMS=M\diag(S);

for t=1:10000
%%    
%%%%%%%   collision     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:lx
    for j=1:ly
        jx=rho(i,j) * u(i,j);jy=rho(i,j) * v(i,j);
        mEq(1,i,j)=  rho(i,j);
        mEq(2,i,j)= -2*rho(i,j) + 3 * (jx^2 + jy^2);
        mEq(3,i,j)=  rho(i,j)   - 3 * (jx^2 + jy^2);
        mEq(4,i,j)=  jx;
        mEq(5,i,j)= -jx;
        mEq(6,i,j)=  jy;
        mEq(7,i,j)= -jy;
        mEq(8,i,j)=  jx^2 - jy^2;
        mEq(9,i,j)=  jx * jy;
    end
end
m = reshape(M * reshape(f,9,lx*ly),9,lx,ly);
f = f - reshape(InvMS * reshape(m-mEq,9,lx*ly),9,lx,ly);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%     streaming      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:9
f(i,:,:) = circshift(f(i,:,:), [0,ex(i),ey(i)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%    Boundary Condition        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Left Wall (bounce back)
f(2,1,:)=f(4,1,:);
f(6,1,:)=f(8,1,:);
f(9,1,:)=f(7,1,:);
%%%% Right Wall (bounce back)
f(4,lx,:)=f(2,lx,:);
f(8,lx,:)=f(6,lx,:);
f(7,lx,:)=f(9,lx,:);
%%%% Bottom Wall (bounce back)
f(3,:,1)=f(5,:,1);
f(6,:,1)=f(8,:,1);
f(7,:,1)=f(9,:,1);
%%%% Top Wall (Zou/He)
rhoLid=f(1,:,ly)+f(2,:,ly)+f(4,:,ly)+2*(f(3,:,ly)+f(7,:,ly)+f(6,:,ly));
f(5,:,ly)=f(3,:,ly);
f(9,:,ly)=f(7,:,ly) + 1/2*(f(4,:,ly)-f(2,:,ly)) + (1/2 * rhoLid * Uo);
f(8,:,ly)=f(6,:,ly) + 1/2*(f(2,:,ly)-f(4,:,ly)) - (1/2 * rhoLid * Uo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%   rho and Velocity    %%%%%%%%%%%%%%%%%%%%%%%%%%%
rho    = reshape(sum(f),lx,ly);
u(:,:) = reshape(ex * reshape(f,9,lx*ly),lx,ly)./rho;
v(:,:) = reshape(ey * reshape(f,9,lx*ly),lx,ly)./rho;
u(1,:) = 0 ; u(lx,:)= 0 ; u(:,1) = 0 ; u(:,ly) = Uo;
v(1,:) = 0 ; v(lx,:)= 0 ; v(:,1) = 0 ; v(:,ly) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(t)
end
%%
%%%%%%%%%%%%%%%   Plotting part   %%%%%%%%%%%%%%%%%%%%%%%%
Psi=zeros(lx,ly);
for i =2:lx-1
    Psi(i,2:ly-1) = Psi(i-1,2:ly-1) - v(i,2:ly-1);
end
contour(0:1/(lx-1):1,0:1/(ly-1):1,Psi',25);
drawnow