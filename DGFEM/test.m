% Simulate the eq: u_t+a u_x=u^2
clear all; close all; %clc;
nx=10;			%discretization in space
p=3;			%polinomial degree
pp=p+1;
stage=1;
rk=stage;		%RK order
T=2;			%Time
cfl=1/(2*p+1);
dx=1/nx;		%Stepwidth in space
dt=dx*cfl;
nt=round(T/dt);
[xl,w]=gauleg(pp);
[P]=legtable(xl,p);
[Pleg]=legtable(xl,p);
x=[];			%for ppotting containing all gauss points of all intervals in space
F_s=zeros(nx,pp,stage);
F_ns=F_s;
tol=[1.d-6,1.d-6]*dx;
 parms = [40,40,-.1,1];
%%%%  Transforming the initial condition to coefficients of Legendre Polinomials  %%%%%%%%%%%%%%%%%%%%%

for i=1:nx
    xi=(2*i-1)*dx/2;      %evaluating the function `func' at the quadrature points
    for m=1:pp
        ffunc(m)=func0(xi+xl(m)*dx/2);
    end
    x=[x,xi+xl*dx/2];
    
    for j=0:p
        u(i,j+1)= sum (ffunc.*P(j+1,:).*w)  * (2*j+1)/2;
    end
end

%% Just for testing the degrees of freedom computation
y0 = (u*P)';
xi = reshape(x,4,nx);
plot(xi,y0);

%%%%%%%%%%  Calculation the Matrix A=int(phi_j + phi_i') , b and c %%%%%%%%
A=zeros(pp,pp);
for i=1:pp
    for j=1:pp
        if j>i && rem(j-i,2)==1
            A(i,j)=2;
        end
    end
    b(i)=(i-1/2)*2/dx;
    c(i)=(-1)^(i-1);
end

%%%%%%%%%%%%%%  RK Method for ODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FS=zeros(nx,pp);
FC=zeros(pp,1);
TIME=0;
TIME = TIME + dt;
%%%%%%%%%%%%  Calculating the d(eta)/d(t) for every timestep i %%%%%%%%%%%%
u_alt = u;
l=1; % stage
% if l==1% Stage 1
% Compute the source term u^2
for i=1:nx
    FC(:)=u_alt(i,:);
    FC=(FC'*Pleg).^2';
    for j=0:p
        FS(i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2;
    end
end
F_ns(1,:,l)=( (A'*u_alt(1,:)' - sum(u_alt(1,1:pp)) + sum(u_alt(nx,1:pp)) * c'+FS(1,1:pp)')' .* b);
for i=2:nx
    F_ns(i,:,l)=( (A'*u_alt(i,:)' - sum(u_alt(i,1:pp)) + sum(u_alt(i-1,1:pp)) * c'+FS(i,1:pp)')' .* b);
end

u = F_ns(:,:,1);
figure
y1 = (u*P)';
xi = reshape(x,4,nx);
subplot(2,1,1); plot(xi,y1,'o-'); title('y1')
subplot(2,1,2); plot(xi,F_ns','o-'); title('F_{ns}')
