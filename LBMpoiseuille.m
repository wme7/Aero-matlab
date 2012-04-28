%*****************Poiseuille flow using LBM***********************%
%Based on procedures as explained in 'Lattice Gas Cellular Automata and Lattice  Boltzmann
%Models'by Wolf Gladrow
%code may have errors as it is my first experience with LBM.Readers are
%suggested to check the code for errors.For feedback
%vinuvargheseijk@gmail.com
OMEGA=0.2;%Relaxation factor
TAU=1/OMEGA;
XMAX=20;%Mesh size in x-direction
YMAX=20;%Mesh size in y-direction
FORCING=(1.024)/YMAX^3;%Forcing term in Poiseuille equation
KVISC = ((1/OMEGA)-0.5)/3;%lattice viscosity
UX=(0.5*FORCING*YMAX*YMAX)/KVISC;%Velocity in x-direction 
cx=[1 0 -1 0 1 -1 -1 1 0];%components of lattice velocities in x-direction
cy=[0 1 0 -1 1 1 -1 -1 0];%components of lattice velocities in y-direction
jx=zeros(XMAX,YMAX);%creating an array for storing momentum in x-direction
jy=zeros(XMAX,YMAX);%creating an array for storing momentum in y-direction
rho=ones(XMAX,YMAX);%creating an array for storing densities
f=zeros(XMAX,YMAX);%creating an array for storing distributions
fprop=zeros(XMAX,YMAX);%creating an array for storing propagating distributions
x=1:XMAX;
y=1:YMAX;
MAXIT=input('Number of iterations');
%Initializataion of distributions
        u=jx./rho;
        v=jy./rho;
      f(x,y,1) = (rho./9)*(1+3.*u+4.5.*u.*u-1.5.*(u.^2+v.^2));
      f(x,y,2) = (rho./9)*(1+3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));
      f(x,y,3) = (rho./9)*(1-3.*u+4.5.*u.*u-1.5*(u.^2+v.^2));
      f(x,y,4) = (rho./9)*(1-3.*v+4.5.*v.*v-1.5*(u.^2+v.^2));

      f(x,y,5) = (rho./36)*(1+3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,6) = (rho./36)*(1+3.*(-u+v)+4.5*(-u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,7) = (rho./36)*(1-3.*(u+v)+4.5*(u+v).^2-1.5*(u.^2+v.^2));
      f(x,y,8) = (rho./36)*(1+3.*(u-v)+4.5*(u-v).^2-1.5*(u.^2+v.^2));

      f(x,y,9) = (4/9)*rho.*(1-1.5*(u.^2+v.^2));
%Initialization completed
%Assign first fprop to equilibrium distributions.
      fprop=f;  
      feq=zeros(20,20,9);
for iter=1:MAXIT
      u=jx./rho;
      v=jy./rho;
                 
      feq(x,2:YMAX-1,1) = (rho(x,2:YMAX-1)./9).*(1+3.*u(x,2:YMAX-1)+4.5.*u(x,2:YMAX-1).*u(x,2:YMAX-1)-1.5.*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,2) = (rho(x,2:YMAX-1)./9).*(1+3.*v(x,2:YMAX-1)+4.5.*v(x,2:YMAX-1).*v(x,2:YMAX-1)-1.5.*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,3) = (rho(x,2:YMAX-1)./9).*(1-3.*u(x,2:YMAX-1)+4.5.*u(x,2:YMAX-1).*u(x,2:YMAX-1)-1.5.*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,4) = (rho(x,2:YMAX-1)./9).*(1-3.*v(x,2:YMAX-1)+4.5.*v(x,2:YMAX-1).*v(x,2:YMAX-1)-1.5.*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));

      feq(x,2:YMAX-1,5) = (rho(x,2:YMAX-1)./36).*(1+3.*(u(x,2:YMAX-1)+v(x,2:YMAX-1))+4.5*(u(x,2:YMAX-1)+v(x,2:YMAX-1)).^2-1.5*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,6) = (rho(x,2:YMAX-1)./36).*(1+3.*(-u(x,2:YMAX-1)+v(x,2:YMAX-1))+4.5*(-u(x,2:YMAX-1)+v(x,2:YMAX-1)).^2-1.5*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,7) = (rho(x,2:YMAX-1)./36).*(1-3.*(u(x,2:YMAX-1)+v(x,2:YMAX-1))+4.5*(u(x,2:YMAX-1)+v(x,2:YMAX-1)).^2-1.5*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,8) = (rho(x,2:YMAX-1)./36).*(1+3.*(u(x,2:YMAX-1)-v(x,2:YMAX-1))+4.5*(u(x,2:YMAX-1)-v(x,2:YMAX-1)).^2-1.5*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      feq(x,2:YMAX-1,9) = ((4/9)*rho(x,2:YMAX-1)).*(1-1.5*(u(x,2:YMAX-1).^2+v(x,2:YMAX-1).^2));
      
%propagating distributions fprop by applying kinetic equation.

      fprop = (1.0-OMEGA).* f + OMEGA.* feq;
      force =FORCING/6;
      fprop(x,2:YMAX-1,1) = fprop(x,2:YMAX-1,1) + force;
      fprop(x,2:YMAX-1,3) = fprop(x,2:YMAX-1,3) - force;
      fprop(x,2:YMAX-1,5) = fprop(x,2:YMAX-1,5) + force;
      fprop(x,2:YMAX-1,6) = fprop(x,2:YMAX-1,6) - force;
      fprop(x,2:YMAX-1,7) = fprop(x,2:YMAX-1,7) - force;
      fprop(x,2:YMAX-1,8) = fprop(x,2:YMAX-1,8) + force;
       
%Copying the boundary values.

            fprop(x,1,:) = f(x,1,:);
            fprop(x,YMAX,:) = f(x,YMAX,:);

%        C6  C2  C5       ^ y
%          \ | /          |
%        C3-C0-C1         | 
%          / | \          | 
%        C7  C4  C8         -----> x
       p=2:XMAX;
       q=2:YMAX;
       r=1:XMAX-1;
       s=1:YMAX-1;
       
%Propagating C1
       
       f(p,y,1) = fprop(p-1,y,1);
      
%Propagating C2

       f(x,q,2) = fprop(x,q-1,2);
      
%Propagating C3

       f(r,y,3) = fprop(r+1,y,3);
      
%Propagating C4
   
       f(x,s,4) = fprop(x,s+1,4);
       
%Propagating C5  

      f(p,q,5) = fprop(p-1,q-1,5);
      
%Propagating C6 

      f(r,p,6) = fprop(r+1,p-1,6);
      
%Propagating C7

      f(r,s,7) = fprop(r+1,s+1,7);
      
%Propagating C8

      f(p,s,8) = fprop(p-1,s+1,8);
      
      
%Propagating C9

      f(x,y,9) = fprop(x,y,9);
      
      
%Complete Bounce Back Boundary Conditions
   
   %1.Implementing Periodic BC.
f(x(1),y,1)=fprop(XMAX,y,1);
f(XMAX,y,3)=fprop(x(1),y,3);
f(x(1),2:YMAX,5)=fprop(XMAX,1:YMAX-1,5);
f(XMAX,2:YMAX,6)=fprop(x(1),1:YMAX-1,6);
f(XMAX,1:YMAX-1,7)=fprop(x(1),2:YMAX,7);
f(x(1),1:YMAX-1,8)=fprop(XMAX,2:YMAX,8);

   %2.Bounce Back Begins.
temp=f(1:XMAX,1,2);
f(1:XMAX,1,2)=f(1:XMAX,1,4);
f(1:XMAX,1,4)=temp;

temp=f(1:XMAX,YMAX,4);
f(1:XMAX,YMAX,4)=f(1:XMAX,YMAX,2);
f(1:XMAX,YMAX,4)=temp;

temp=f(1:XMAX,1,5);
f(1:XMAX,1,5)=f(1:XMAX,1,7);
f(1:XMAX,1,7)=temp;

temp=f(1:XMAX,YMAX,7);
f(1:XMAX,YMAX,7)=f(1:XMAX,YMAX,5);
f(1:XMAX,YMAX,5)=temp;

temp=f(1:XMAX,1,6);
f(1:XMAX,1,6)=f(1:XMAX,1,8);
f(1:XMAX,1,8)=temp;

temp=f(1:XMAX,YMAX,8);
f(1:XMAX,YMAX,8)=f(1:XMAX,YMAX,6);
f(1:XMAX,YMAX,6)=temp;

rho(x,y) = f(x,y,1)+f(x,y,2)+f(x,y,3)+f(x,y,4)+f(x,y,5)+f(x,y,6)+f(x,y,7)+f(x,y,8)+f(x,y,9);
jx(x,y)=f(x,y,1)-f(x,y,3)+f(x,y,5)-f(x,y,6)-f(x,y,7)+f(x,y,8);%Distributions multiplied by lattice velocities in x-directions
jy(x,y)=f(x,y,2)-f(x,y,4)+f(x,y,5)+f(x,y,6)-f(x,y,7)-f(x,y,8);%Distributions multiplied by lattice velocities in y-directions
uprofile(1:YMAX)=sum(jx(1:XMAX,1:YMAX)./rho(1:XMAX,1:YMAX))/XMAX;%Velocity profile in x-direction
plot(uprofile)
pause(.1)
end
