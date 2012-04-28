%GLBE based incompressible MRT model for gravity 
%driven fluid flow simulation in a parallel channel.
%
%D2Q9: 2D- 9 velocity model

%lx is length of domain and ly is width. 
%g is applied gravity in X-direction
%Copyright: Shadab Anwar, Florida International University, USA.
%sanwa001@fiu.edu
%*********************************************************************
clear all;close all;
wa(9)=4/9;wa(1:4)=1/9;wa(5:8)=1/36;
lx=8;ly=33; rho0_in=1.0001; rho0_out=1.;
u=zeros(ly,lx,2); g= 1e-04;
u_pois=zeros(1,ly); f=zeros(ly,lx,9);
ftemp=zeros(ly,lx,9);rho=ones(ly,lx);
u(:,:,1)=0.0041;
M=[   1,1,1,1,1,1,1,1,1;
     -1,-1,-1,-1,2,2,2,2,-4;
     -2,-2,-2,-2,1,1,1,1,4;
      1,0,-1,0,1,-1,-1,1,0;
     -2,0,2,0,1,-1,-1,1,0;
      0,1,0,-1,1,1,-1,-1,0;
      0,-2,0,2,1,1,-1,-1,0;
      1,-1,1,-1,0,0,0,0,0;
      0,0,0,0,1,-1,1,-1,0;
];

s=[0,-1.64,-1.54,0,-1.9,0,-1.9,-1.,-1.]; % relaxation parameters
%s(9) controls the kinematic viscosity, and hence the Reynolds number.

ts=300;

a2=-8;a3=4;c1=-2;g1=2/3;g3=2/3;g2=18;g4=-18;  % a3 and g4 are free to tune
e(1,:)=[1,0,-1,0,1,-1,-1,1,0]; %ex
e(2,:)=[0,1,0,-1,1,1,-1,-1,0]; %ey

%i=10;j=5;
%rho=1;

 tau(1)=-1/s(9); 
 nu=(1/3)*(tau(1)-0.5);
 ni=ly;
 %dpdl=(1/3)*(1e-04)/lx;
  mm=1;

for i=-(ni-1)/2:(ni-2)/2
u_pois(mm)=(g/(2*nu))*(((ni-2)/2)^2-i^2);
%u_pois(mm)=(dpdl/(2*nu))*(((ni-2)/2)^2-i^2);
u_pois(1)=0.0;
mm=mm+1;
end
u_pois(mm)=0.0;

%initialize f's
for i=1:lx
    for j=1:ly

u(j,i,1)=u_pois(j);

jx=u(j,i,1);
jy=u(j,i,2);

m_eq(1) = rho(j,i); %density
m_eq(2) = (a2/4)*rho(j,i) + (1/6)*g2*(jx^2+jy^2);  %energy
m_eq(3) = (a3/4)*rho(j,i) + (1/6)*g4*(jx^2+jy^2);  %energy square
m_eq(4) = jx; %momentum in x-dir
m_eq(5) = (1/2)*c1*jx; %energy flux in x-dir
m_eq(6) = jy; %momentum in y-dir
m_eq(7) = (1/2)*c1*jy; %energy flux in y-dir
m_eq(8) = (3/2)*g1*(jx^2 - jy^2); %diagonal comp of stress tensor
m_eq(9) = (3/2)*g3*jx*jy;  %off diagonal comp of stress tensor
            
f(j,i,:)=inv(M)*(m_eq');
    end
end

        F(9) = 3*wa(9)*( e(1,9)*g );
        for a=1:8
             F(a)=3*wa(a)*( e(1,a)*g );  %F is applied body force (gravity)
        end
for t=1:ts
    t
%Macroscopic variable
for i=1:lx
    for j=1:ly
        
        rho(j,i)=sum(f(j,i,:));
        u(j,i,:)=0.;
        
         for a=1:9
        u(j,i,1) = u(j,i,1) + e(1,a)*f(j,i,a);
        u(j,i,2) = u(j,i,2) + e(2,a)*f(j,i,a);
          end
          
      u(j,i,1)= u(j,i,1)/rho(j,i);
      u(j,i,2)= u(j,i,2)/rho(j,i);
  
f_pre=ones(9,1);f_post=ones(9,1);

%Compute Meq or moment of feq
jx=u(j,i,1);
jy=u(j,i,2);

m_eq(1) = rho(j,i); %density
m_eq(2) = (a2/4)*rho(j,i) + (1/6)*g2*(jx^2+jy^2);  %energy
m_eq(3) = (a3/4)*rho(j,i) + (1/6)*g4*(jx^2+jy^2);  %energy square
m_eq(4) = rho(j,i)*jx; %momentum in x-dir
m_eq(5) = (1/2)*c1*jx; %energy flux in x-dir
m_eq(6) = rho(j,i)*jy; %momentum in y-dir
m_eq(7) = (1/2)*c1*jy; %energy flux in y-dir
m_eq(8) = (3/2)*g1*(jx^2 - jy^2); %diagonal comp of stress tensor
m_eq(9) = (3/2)*g3*jx*jy;  %off diagonal comp of stress tensor

%Collision operator
            if (j==1 || j==ly) % standard BOUNCE-BACK on Solid nodes
    
    temp = f(j,i,1); f(j,i,1) = f(j,i,3); f(j,i,3) = temp;
    temp = f(j,i,2); f(j,i,2) = f(j,i,4); f(j,i,4) = temp;
    temp = f(j,i,5); f(j,i,5) = f(j,i,7); f(j,i,7) = temp;
    temp = f(j,i,6); f(j,i,6) = f(j,i,8); f(j,i,8) = temp;
            else
f_pre(1:9)=f(j,i,:);

m=M*f_pre;
dm=m+diag(s)*(m-m_eq');%+diag(ss)*m_eq';
df=inv(M)*dm;

f(j,i,:)=df'+F;
            end %bounceback

   end %j loop
end %i loop


%streaming , PERIODIC BOUNDARY ALONG X-AXIS
for   i=1:lx
          if(i>1)
              in=i-1;
          else
              in=lx;
          end
    
          if(i<lx)
              ip=i+1;
          else
              ip=1;
          end
          
        for j=1:ly
      if(j>1)
      jn=j-1;
      else
      jn=ly;
      end

      if(j<ly)
      jp=j+1;
      else
      jp=1;
      end

    ftemp(j,i,9)   =  f(j,i,9);
    ftemp(j,ip,1)  =  f(j,i,1);
    ftemp(jp,i,2)  =  f(j,i,2);
    ftemp(j,in,3)  =  f(j,i,3);
    ftemp(jn,i,4)  =  f(j,i,4);
    ftemp(jp,ip,5) =  f(j,i,5);
    ftemp(jp,in,6) =  f(j,i,6);
    ftemp(jn,in,7) =  f(j,i,7);
    ftemp(jn,ip,8) =  f(j,i,8);
    
        end % streaming j loop
     
 end % streaming i loop
 
    f=ftemp; 

end % time loop

ux=u(:,:,1);
uy=u(:,:,2);


figure
plot(uy(1:ly,lx/2))
title('Velocity profile')
%figure
hold on; 
 plot(u_pois(1:ly),'r')
 
 plot(ux(1:ly,lx/2),'o')
 legend('LBM(uy)','Poiseuille (ux)','LBM(ux)',0)
 %title('ux')
 hold off
 

 
