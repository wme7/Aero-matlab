%------------------------------------------------------------------------%
% Lax-Friedrichs method to solve 1-D Euler equations
% -----------------------------------------------------------------------%
clear;
clc;
close all;
n=601;              %Number of grid points
L=10;               %Length of domain
h=L/(n-1);          %Spatial step size
CFL=0.2;            %CFL number for stability
alpha=0.1           %Parameter for artificial viscosity
t_final=3.9e-3;     %Final time
x=0:h:L;
gamma=1.4;          %Ratio of specific heats for ideal di-atomic gas
%------------------------------------------------------------------------%
% Define initial conditions
%------------------------------------------------------------------------%
p_l=1e5;            %Pressure in left side of shock tube at t=0
p_r=1e3;            %Pressure in right side of shock tuve at t=0
rho_l=1;            %Density at left side of shock tube at t=0
rho_r=0.01;         %Density at right side of shock tube at t=0
u_l=0;              %Velocity in left side of shock tube at t=0
u_r=0;              %Velocity in right side of shock tube at t=0
p(1:1:(n+1)/2)=p_l;
p((n+3)/2:1:n)=p_r;
rho(1:1:(n+1)/2)=rho_l;
rho((n+3)/2:1:n)=rho_r;
u(1:1:(n+1)/2)=u_l;
u((n+3)/2:1:n)=u_r;
E=p./((gamma-1)*rho)+0.5*u.^2;  %Total Energy
a=sqrt(gamma*p./rho);           %Speed of sound
dt=CFL*h/max(abs(a));           %Time step based on CFL number
step=0;
%------------------------------------------------------------------------%
% Time integration begins
%------------------------------------------------------------------------%
for t=dt:dt:t_final
    %Define q & F matrix
    q=[rho; rho.*u; rho.*E];
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)];
    %Calculate q* and flow parameters
    q_star(1:3,1:n-1)=0.5*(q(1:3,1:n-1)+q(1:3,2:n))-dt/(2*h)*(F(1:3,2:n)-F(1:3,1:n-1));
    rho(1:n-1)=q_star(1,1:n-1);
    u(1:n-1)=q_star(2,1:n-1)./rho(1:n-1);
    E(1:n-1)=q_star(3,1:n-1)./rho(1:n-1);
    p(1:n-1)=(gamma-1)*rho(1:n-1).*(E(1:n-1)-0.5*u(1:n-1).^2);
    %Calculate F*
    F_star(1:3,1:n-1)=[rho(1:n-1).*u(1:n-1); rho(1:n-1).*u(1:n-1).^2+p(1:n-1); u(1:n-1).*(rho(1:n-1).*E(1:n-1)+p(1:n-1))];
    mx=[zeros(1,n-1);ones(1,n-1);u(1:n-1)];
    %Calculate arfiticial viscosity
    for i=1:3
        visc(i,1:n-1)=alpha*h^2*rho(1:n-1).*abs((u(2:n)-u(1:n-1))/h).*((u(2:n)-u(1:n-1))/h).*mx(i,1:n-1);
    end
    %Update F* and q matrix
    F_star=F_star-visc;
    q(1:3,2:n-1)=q(1:3,2:n-1)-dt/h*(F_star(1:3,2:n-1)-F_star(1:3,1:n-2));
    rho=q(1,1:n);
    u=q(2,1:n)./rho(1:n);
    E=q(3,1:n)./rho;
    p=(gamma-1)*rho.*(E-0.5*u.^2);
    step=step+1;
end
%calculation of flow parameters
a=sqrt(gamma*p./rho);
M=u./a;
p_ref=101325;           %Reference air pressure (N/m^2)
rho_ref=1.225;          %Reference air density (kg/m^3)
s=1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho)); %Entropy w.r.t reference condition
Q=rho.*u;               %Mass Flow rate per unit area
%------------------------------------------------------------------------%
offset=0.05;
subplot(231);plot(x,p,'k');xlabel('X-Coordinate (m)');ylabel('Pressure (Pa)');ylim([min(p)-(offset)*max(p) (1+offset)*max(p)]);
subplot(232);plot(x,s,'k');xlabel('X-Coordinate (m)');ylabel('Entropy/R gas');ylim([min(s)-(offset)*max(s) (1+offset)*max(s)]);
subplot(233);plot(x,u,'k');xlabel('X-Coordinate (m)');ylabel('Velocity (m/s)');ylim([min(u)-(offset)*max(u) (1+offset)*max(u)]);
subplot(234);plot(x,M,'k');xlabel('X-Coordinate (m)');ylabel('Mach number');ylim([min(M)-(offset)*max(M) (1+offset)*max(M)]);
subplot(235);plot(x,rho,'k');xlabel('X-Coordinate (m)');ylabel('Density (kg/m^3)');ylim([min(rho)-(offset)*max(rho) (1+offset)*max(rho)]);
subplot(236);plot(x,Q,'k');xlabel('X-Coordinate (m)');ylabel('Mass Flow (kg/m^2s)');ylim([min(Q)-(offset)*max(Q) (1+offset)*max(Q)]);
%------------------------------------------------------------------------%
%END ;D

%%
% 
% $$e^{\pi i} + 1 = 0$$
% 
