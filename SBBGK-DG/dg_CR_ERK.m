
% Simulate the eq: u_t+a u_x=u^2
% function [error,error2]=dg_transport_RK(nx,nt,p,rk)
clear all

coeffi_RK % Calculation of the coefficients alpha for ARK
gamma=const_a_I(2,1);

nx=160;			%discretization in space
p=3;			%polinomial degree
pp=p+1;
stage=6;
rk=stage;			%RK order
T=0.25;			%Time
cfl=1/(2*p+1);

%dt=T/nt;		%stepwidth in time
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
%%%%%%%%%%%%%%  Transforming the initial condition to coefficients of Legendre Polinomials  %%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%  Calculation the Matrix A=int(phi_j + phi_i') , b and c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(pp,pp);
for i=1:pp
    for j=1:pp
        if j>i& rem(j-i,2)==1
            A(i,j)=2;
        end
    end
    b(i)=(i-1/2)*2/dx;
    c(i)=(-1)^(i-1);
end


%%%%%%%%%%%%%%  Calculation of the coefficients alpha for SSP-RK  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% alpha(1)=1;
% 
% for m=1:rk
%     for k=(m-1):(-1):1
%         alpha(k+1)=1/k*alpha(k);
%     end
%     alpha(m)=1/factorial(m);
%     alpha(1)=1-sum(alpha(2:m));
% end


%%%%%%%%%%%%%%  RK Method for ODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FS=zeros(nx,pp);
FC=zeros(pp,1);
TIME=0;
ISTOP=0;
uold=u;
while  ISTOP ==0   
    
    if (TIME +dt> T)
        dt=T- TIME;
        ISTOP = 1;
    elseif TIME+dt==T
        ISTOP = 1;
    else
        ISTOP = 0;
    end
    TIME = TIME + dt;
    u_alt=u;    
    for l=1:rk        
        %%%%%%%%%%%%  Calculating the d(eta)/d(t) for every timestep i %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if l==1% Stage 1
            
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
                F_ns(i,:,l)=( ( A'*u_alt(i,:)' - sum(u_alt(i,1:pp)) + sum(u_alt(i-1,1:pp)) * c'+FS(i,1:pp)')' .* b);
            end
        else
             % Stages 2-6          
            
            for i=1:nx
                % Compute the source term u^2
                FC(:)=u_alt(i,:);
                FC=(FC'*Pleg).^2';
                for j=0:p
                    FS(i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2;
                end
            end
            % Right-going wave u_t+u_x=u^2
            
            F_ns(1,:,l)=( (A'*u_alt(1,:)' - sum(u_alt(1,1:pp)) + sum(u_alt(nx,1:pp)) * c'+FS(1,1:pp)')' .* b);
            %  terms related to sum(u_alt) are related the numerical flux
            %  (upwinding flux)
            for i=2:nx
                F_ns(i,:,l)=( ( A'*u_alt(i,:)' - sum(u_alt(i,1:pp)) + sum(u_alt(i-1,1:pp)) * c'+FS(i,1:pp)')' .* b);
            end
        end      
        
        if l<stage
            u = u + dt*const_b(l)*(F_s(:,:,l)+F_ns(:,:,l));
            u_alt = uold;
            for j=1:l %u_alt=Un+Xi
                u_alt = u_alt + dt*(const_a_I(l+1,j)*F_s(:,:,j) + const_a_E(l+1,j)*F_ns(:,:,j));
            end
        else
            % Final Stage and update the new U, Eq. 42
            u = u + dt*const_b(l)*(F_s(:,:,l)+F_ns(:,:,l));
        end        
    end
    uold=u;
    y=u*P;
    z=[];
    w_collapse=[];
    for i=1:size(y,1)
        z=[z,y(i,:)];
        w_collapse=[w_collapse,w];
    end
    plot(x,z);
    drawnow
end


%%%%%%%%%%%%%  plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error=sqrt( sum(((1./(1./func0(x-T)-T)-z).^2).*w_collapse)*dx/2 )    % Error in L2 Norm
error2=max(abs((1./(1./func0(x-T)-T)-z)))				   % Error in max Norm


% error=sqrt( sum(((func0(x-T).*exp(T)-z).^2).*w_collapse)*dx/2 )    % Error in L2 Norm
% error2=max(abs((func0(x-T).*exp(T)-z)))				   % Error in max Norm
