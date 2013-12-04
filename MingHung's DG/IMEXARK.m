% PDETIMEDEMO
% This program creates the left preconditionied time-dependent PDE example
% from Chapter 3
%
%
% The example has the form
%
% u_t = - L(u) + f, 0 < x,y < 1
% with homogeneous Dirichlet BC
%
% In terms of the program below, the equation is
% Lu = - (u_xx + u_yy)
% f= 4pi^2 exp(-4 pi^2) sin(2pi x) sin (2pi y)
% u= exp(-4 pi^2) sin(2pi x) sin (2pi y)
% Within the nonlinear map, pdetime.m, we apply fish2d.m as
% a preconditioner. 
%
% The action of the Laplacian, partial wrt x, and partial wrt y
% are computed in the routines lapmf, dxmf, and dymf. Here mf
% stands for MATRIX FREE.
%
global rhsf uold dt u_alt gamma;

coeffi_RK
gamma=const_a_I(2,1);
stage=6;
m=64*2-1; 
n=m*m;
h=1/(m+1);
h2=(m+1)*(m+1);
h1=(m+1);
tolh=1/h2;
mres=3;
%
% set up the equation
%
x=1:m;
e=x';
x=x'*h;
%
% Form the steady-state solution: sol(x,y)= 10 x y (1-x) (1-y) exp(x^4.5)
%
% zsoly=x.*(1-x);
% zsolx=zsoly.*exp(x.^4.5);
% solt=10*zsolx*zsoly';
% sol=solt(:);
zsoly=sin(2*pi*x);
zsolx=sin(2*pi*x);
solt=zsolx*zsoly';
sol=solt(:);

clear rhs; clear axt; clear at; clear ayt; 
%
% Fix the right side so that the steady-state solution is known
%
%rhsf = lapmf(sol) + 20*sol.*(dxmf(sol)+dymf(sol));
rhsf =0;% 4*pi^2*exp(-4*pi^2)*sin(2*pi*x)*sin*(2*pi*y);
%
% Set the initial data to zero. 
%
u0=zeros(m*m,1);
u0=solt(:);
uold=u0;
uinit=uold;

nt=8;
dt=0.1/nt;

F_s  = zeros(length(u0),6);

%
%
% Use nsoli.m with GMRES as the solver. 
% Converged result from previous time step is initial iterate.
%
tol=[1.d-6,1.d-6]*tolh;
x=zeros(n,1); parms = [40,40,-.1,1];
t=0;
for it=1:nt
    
    sol=exp(-4*pi^2*t)*solt(:);
    norm(uold-sol);
    u_num=uold;
    u_alt = uold;
    %u_n = uold;
    for i=1:stage
        % The equation is u_t = (u_xx + u_yy) - 20*u*(u_x + u_y)  + f
    if i==1
        %F_s(:,i)  = -lapmf(u_alt)+4*pi^2*u_alt;
        F_ns(:,i)  = 4*pi^2*u_alt;    
        F_s(:,i)  = -lapmf(u_alt);     
        %F_s(:,i)  = -lapmf(u_alt)-20*u_alt.*(dxmf(u_alt)+dymf(u_alt))+rhsf;
%        F_s(:,i)  = Luq*v_alt;
%         F_ns(:,i) = -current(v_alt,m_alt,H_alt,J_alt,d_alt,f_alt,...
%                          X_alt,Ca_alt) + Stim;
    else
        [u_n, it_histg, ierr] = nsoli(uold,'pdetimeimex',tol,parms);
         %F_s(:,i)  = -lapmf(u_n)+4*pi^2*u_n;
         F_ns(:,i)  =   4*pi^2*u_n;
         F_s(:,i)  = -lapmf(u_n);
 %        F_s(:,i)  = -lapmf(u_alt)-20*u_alt.*(dxmf(u_alt)+dymf(u_alt))+rhsf;
%        F_s(:,i)  = Luq*inv_mat_IQ*v_alt;
%         F_ns(:,i) = -current(inv_mat_IQ*v_alt,m_alt,H_alt,J_alt,d_alt,f_alt,...
%                          X_alt,Ca_alt) + Stim;
    end    
    if i<stage
        %u_num = u_num + dt*const_b(i)*(F_s(:,i)); %+F_ns(:,i));      
        u_num = u_num + dt*const_b(i)*(F_s(:,i)+F_ns(:,i));    
        u_alt = uold;
%         v_num = v_num + dt*const_b(i)*(F_s(:,i)+F_ns(:,i));        
%         v_alt = v_n;
        for j=1:i %u_alt=Un+Xi
            %u_alt = u_alt + dt*(const_a_I(i+1,j)*F_s(:,j) );%+ const_a_E(i+1,j)*F_ns(:,j));  
            u_alt = u_alt + dt*(const_a_I(i+1,j)*F_s(:,j) + const_a_E(i+1,j)*F_ns(:,j));  
%            v_alt = v_alt + dt*(const_a_I(i+1,j)*F_s(:,j) + const_a_E(i+1,j)*F_ns(:,j));  
        end
    else
        %u_num = u_num + dt*const_b(i)*(F_s(:,i));%+F_ns(:,i)); 
        u_num = u_num + dt*const_b(i)*(F_s(:,i)+F_ns(:,i));
%        v_num = v_num + dt*const_b(i)*(F_s(:,i)+F_ns(:,i));
    end
end
    
    %[unew, it_histg, ierr] = nsoli(uinit,'pdetime',tol,parms);

    unew=u_num;
    uinit=uold;
    uold=unew;
    t=t+dt;
end
sol=exp(-4*pi^2*t)*solt(:);
norm(uold-sol)
max(abs(uold-sol))
