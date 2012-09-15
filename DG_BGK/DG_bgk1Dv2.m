%% Program Direct solver using DG and BGK approximation
% Matlab Modification of lab007's Direct Solver Code

% clear memory and figures
clear all; close all; clc;

% start counting time
tic 

% Parameters
nx = 10; % number of elements
p  = 5;	 % polinomial degree
pp = p+1;
rk = pp; % RK order
CFL= 1/(2*p+1);
bb = 1;

GHNC  = 0;
%CFL  = 0.9;
OUTTIME = 0.1;
TAU	  = 0.01;% !RELAXATION TIME

IT = 1;

NV = 20;
NVh=20/2;
GH =[-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738, ...
    -2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,...
    -0.737473728545,-0.245340708301,0.245340708301,0.737473728545,...
    1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,...
    3.94476404012,4.60368244955,5.38748089001];
wp  =[0.898591961453,0.704332961176,0.62227869619,0.575262442852,...
    0.544851742366,0.524080350949,0.509679027117,0.499920871336,...
    0.493843385272,0.490921500667,0.490921500667,0.493843385272,...
    0.499920871336,0.509679027117,0.524080350949,0.544851742366,...
    0.575262442852,0.62227869619,0.704332961176,0.898591961453];

V=-GH;

%% Initial Conditions (IC):
UL  = 0.;
TL  = 4.38385;
ZL  = 0.2253353;
UR  = 0.;
TR  = 8.972544;
ZR  = 0.1204582;
% UR  = UL;
% TR  = TL;
% ZR  = ZL;

dx=1/nx;    %Stepwidth in space
dt=dx*CFL/60; %abs(V(1))

nt=round(OUTTIME/dt);
[xl,w]=gauleg(pp);
[Pleg]=legtable(xl,p);

%% Initialize variables;
F=zeros(NV,nx,pp);
FEQ=zeros(NV,nx,pp);
F_tmp=zeros(NV,nx,pp);
F_new=zeros(NV,nx,pp);
FS=zeros(NV,nx,pp);
F_loc=zeros(NV,pp);
SR=zeros(nx,pp);
SU=zeros(nx,pp);
SE=zeros(nx,pp);
SAV=zeros(nx,pp);
R=zeros(nx,pp);
P=zeros(nx,pp);
U=zeros(nx,pp);
T=zeros(nx,pp);
Z=zeros(nx,pp);
ET=zeros(nx,pp);
AV=zeros(nx,pp);
x=zeros(1,nx*pp);
ffunc=zeros(1,pp);
alpha=zeros(1,rk);

%% Transforming the IC to coefficients of Legendre Polinomials

for i=1:nx
    xi=(2*i-1)*dx/2; %evaluating the function `func' at the quadrature points
    x((i-1)*pp+1:i*pp)=xi+xl*dx/2;
    for m=1:pp
        if(xi+xl(m)*dx/2 <= 0.5)
            U(i,m) = UL;
            T(i,m) = TL;
            Z(i,m) = ZL;
        else
            U(i,m) = UR;
            T(i,m) = TR;
            Z(i,m) = ZR;
        end
        
    end
    for K = 1: NV
        for m=1:pp
            ffunc(m)  = 1/((exp((V(K)-U(i,m))*(V(K)-U(i,m))/T(i,m))/Z(i,m)) + IT);
            %ffunc(m)=func(xi+xl(m)*dx/2);
        end
        for j=0:p
            F(K,i,j+1)= sum (ffunc.*Pleg(j+1,:).*w)*(2*j+1)/2;
        end
    end
end

%%
for i=1:nx
    Mtemp=zeros(NV,pp);
    for K=1:NV
        Mtemp(K,:)=F(K,i,:);
    end
    F_loc(:,:)=Mtemp*Pleg;
    
    for m=1:pp
        SR(i,:) = wp * F_loc;
        SU(i,m) = sum(wp.*F_loc(:,m)'.* V);
        SE(i,m) = sum(wp.*F_loc(:,m)'.* V.^2)/2;
        SAV(i,m)= sum(wp.*F_loc(:,m)'.* abs(V));
    end
end
R    = SR;
U   = SU./SR;
ET  = SE;
AV  = SAV;

% fp_plot=zeros(nx*pp,1);
% fm_plot=zeros(nx*pp,1);
% Kplot=1;
% for i=1:nx
%     fp_plot((i-1)*pp+1:i*pp)=Floc(Kplot,i,1:pp);
%     fm_plot((i-1)*pp+1:i*pp)=Floc(NV-Kplot+1,i,1:pp);
% end
% %

r_plot=reshape(R',nx*pp,1);
u_plot=reshape(U',nx*pp,1);
if bb==1
    %     subplot(121),
    wave_handleu=plot(x,r_plot);
    axis([-0.2, 1.2, 0.5, 0.85]);
    xlabel('x'); ylabel('R(x,t)')
    set(wave_handleu,'EraseMode','Xor'); drawnow
    %     subplot(122),wave_handlev=plot(x,u_plot);
    %     axis([-0.5, 1.5, -1, 1]);
    %     xlabel('x'); ylabel('U(x,t)')
    %     set(wave_handlev,'EraseMode','Xor'); drawnow
end

% if bb==1
%     figure(1)
%     title(sprintf('Iteration %d',ITER));
%     subplot(211),wave_handleu=plot(x,fp_plot);
%     axis([0, 1, 0, 2]);
%     xlabel('x'); ylabel('u(x,t)')
%     set(wave_handleu,'EraseMode','Xor'); drawnow
%     subplot(212),wave_handlev=plot(x,fm_plot);
%     axis([0, 1, 0, 2]);
%     xlabel('x'); ylabel('v(x,t)')
%     set(wave_handlev,'EraseMode','Xor'); drawnow
% end

%%
pause(1)

A=zeros(pp,pp);
b=zeros(1,pp);
c=zeros(1,pp);
for i=1:pp
    for j=1:pp
        if j>i && rem(j-i,2)==1
            A(i,j)=2;
        end
    end
    b(i)=(i-1/2)*2/dx;
    c(i)=(-1)^(i-1);
end
alpha(1)=1;
for m=1:rk
    for k=(m-1):(-1):1
        alpha(k+1)=1/k*alpha(k);
    end
    alpha(m)=1/factorial(m);
    alpha(1)=1-sum(alpha(2:m));
end

ITER  = 1;
TIME  = 0;
ISTOP = 0;
FC=zeros(pp,1);
FN=zeros(pp,1);
FB=zeros(pp,1);


while ISTOP ==0
    VIS(1:nx,1:pp) = TAU;
    
    %dt = dx * CFL/V(1);
    TIME = TIME + dt;
    dtdx = dt/dx;
    
    if (TIME > OUTTIME)
        DTCFL = OUTTIME - (TIME - dt) ;
        TIME = OUTTIME;
        dt = DTCFL;
        dtdx = dt /dx;
        ISTOP = 1;
    end
    
    %         F(K,i,j+1)=phi;
    %         phi=alpha(1)*phi;
    %         psi_alt=psi;
    %         psi=alpha(1)*psi;
    
%% Calculating the d(eta)/d(t) for every timestep i 
    F_tmp=F;
    F_new=alpha(1)*F;
    
    for l=1:rk
        
        for i = 1:nx
            for K = 1:NV
                for m=1:pp
                    FEQ(K,i,m)   = 1/((exp( (V(K)-U(i,m))^2 /T(i,m))/Z(i,m)) + IT );
                end
            end
        end
        
        for i=1:nx
            for K = 1: NV
                FC(:)=F(K,i,:);
                FB(:)=FEQ(K,i,:);
                FC=(FC'*Pleg-FB')';
                
                for j=0:p
                    FS(K,i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2/VIS(i,j+1);
                end
            end
            
            
        end
        for K=1:NVh
            %phi_t(1,:)=( (A'*phi_alt(1,:)' - sum(phi_alt(1,1:pp))- sum(psi_alt(1,1:pp).* c) * c')' .* b); %BC reflecting
            FC(:)=F(K,1,:);
            FB(:)=F(NV-K+1,1,:);
            FN(:)=FS(K,1,:);
            F_tmp(K,1,:)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FB'.* c) * c'-FN)' .* b);
            for i=2:nx
                FN(:)=F(K,i-1,:);
                FC(:)=F(K,i,:);
                FB(:)=FS(K,i,:);
                F_tmp(K,i,:)=( (V(K)*A'*FC -V(K)* sum(FC) +V(K)* sum(FN) * c'-FB)' .* b);
                % phi_t(i,:)=( (A'*phi_alt(i,:)' - sum(phi_alt(i,1:pp)) + sum(phi_alt(i-1,1:pp)) * c')' .* b);
            end
            
            for i=1:nx-1
                FN(:)=F(NVh+K,i+1,:);
                FC(:)=F(NVh+K,i,:);
                FB(:)=FS(NVh+K,i,:);
                F_tmp(NVh+K,i,:)=( (-V(NVh+K)*A'*(-FC) -V(NVh+K)* sum(FN'.* c) + V(NVh+K)*sum(FC'.* c) * c'-FB)' .* b);
                %F_tmp(NVh+K,i,:)=( (A'*-psi_alt(i,:)' + sum(psi_alt(i+1,1:pp).* c) - sum(psi_alt(i,1:pp).* c) * c')' .* b);
            end
            %            psi_t(nx,:)=( (A'*-psi_alt(nx,:)' - sum(phi_alt(nx,1:pp)) - sum(psi_alt(nx,1:pp).*c) * c')' .* b);
            FC(:)=F(NVh+K,nx,:);
            FB(:)=F(NVh-K+1,nx,:);
            FN(:)=FS(NVh+K,nx,:);
            F_tmp(NVh+K,nx,:)=( (-V(NVh+K)*A'*(-FC) -V(NVh+K)* sum(FB) +V(NVh+K)* sum(FC'.*c) * c'-FN)' .* b);
            
        end % loop for NV
       
        if l<rk
            %                 phi=phi+ alpha(l+1)*(phi_alt+dt*phi_t);
            %                 phi_alt=phi_alt+dt*phi_t;
            %                 psi=psi+ alpha(l+1)*(psi_alt+dt*psi_t);
            %                 psi_alt=psi_alt+dt*psi_t;
            F_new=F_new+ alpha(l+1)*(F+dt*F_tmp);
            F=F+dt*F_tmp;
            
            for i=1:nx
                Mtemp=zeros(NV,pp);
                for K=1:NV
                    Mtemp(K,:)=F(K,i,:);
                end
                F_loc(:,:)=Mtemp*Pleg;
                SR(i,:) = wp * F_loc;
                for m=1:pp
                    SU(i,m) = sum(wp.*F_loc(:,m)'.* V);
                    SE(i,m) = sum(wp.*F_loc(:,m)'.* V.^2)/2;
                    SAV(i,m)= sum(wp.*F_loc(:,m)'.* abs(V));
                    %         for K = 1: NV
                    %             SR(i,m) = SR(i,m) + wp(K) * Floc(K,i,m);
                    %             SU(i,m) = SU(i,m) + wp(K) *  Floc(K,i,m) * V(K);
                    %             SE(i,m) = SE(i,m) + wp(K) *  Floc(K,i,m) * (0.5 * V(K) * V(K));
                    %             SAV(i,m) = SAV(i,m) + wp(K) *  Floc(K,i,m) * abs(V(K));
                    %         end
                    %                     R(i,m)    = SR(i,m);
                    %                     U(i,m)    = SU(i,m)/SR(i,m);
                    %                     ET(i,m)   = SE(i,m);
                    %                     AV(i,m)   = SAV(i,m);
                end
                
            end
            R    = SR;
            U   = SU./SR;
            ET  = SE;
            AV  = SAV;
            switch IT
                case 0
                    %     if (IT == 0)
                    for i=1:nx
                        for m=1:pp
                            T(i,m)    = 4*ET(i,m)/R(i,m) - 2*U(i,m)^2;
                            Z(i,m)    = R(i,m) / SQRT(pi* T(i,m));
                            P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                        end
                    end
                case -1 %(IT == -1)
                    for i=1:nx
                        for m=1:pp
                            ZA = 0.0001;
                            ZB = 0.99;
                            GA12 = 0;
                            GB12 = 0;
                            GA32 = 0;
                            GB32 = 0;
                            for L = 1:50
                                GA12 = GA12 + (ZA^L)/(L^0.5);
                                GB12 = GB12 + (ZB^L)/(L^0.5);
                                GA32 = GA32 + (ZA^L)/(L^1.5);
                                GB32 = GB32 + (ZB^L)/(L^1.5);
                            end
                            PSIA = 2*ET(i,m) - GA32*(R(i,m)/GA12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                            PSIB = 2*ET(i,m) - GB32*(R(i,m)/GB12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                            while (abs(ZA-ZB) > 0.00001)
                                
                                ZC = (ZA + ZB)/2;
                                GC12 = 0;
                                GC32 = 0;
                                GC52 = 0;
                                for L = 1:50
                                    GC12 = GC12 + (ZC^L)/(L^0.5);
                                    GC32 = GC32 + (ZC^L)/(L^1.5);
                                    GC52 = GC52 + (ZC^L)/(L^2.5);
                                end
                                PSIC = 2*ET(i,m) - GC32*(R(i,m)/GC12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                                
                                if ((PSIA*PSIC) < 0)
                                    ZB = ZC;
                                    GB12 = GC12;
                                    GB32 = GC32;
                                    PSIB =PSIC;
                                else
                                    ZA = ZC;
                                    GA12 = GC12;
                                    GA32 = GC32;
                                    PSIA =PSIC;
                                end
                            end
                            Z(i,m) = ZC;
                            T(i,m) = R(i,m)^2 / (pi*GC12^2 );
                            P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                        end
                    end
                case 1            %(IT == 1)
                    for i=1:nx
                        for m=1:pp
                            ZA = 0.0001;
                            ZB = 0.99;
                            GA12 = 0;
                            GB12 = 0;
                            GA32 = 0;
                            GB32 = 0;
                            for L = 1:50
                                GA12 = GA12 + (ZA^L)*(-1)^(L-1)/(L^0.5);
                                GB12 = GB12 + (ZB^L)*(-1)^(L-1)/(L^0.5);
                                GA32 = GA32 + (ZA^L)*(-1)^(L-1)/(L^1.5);
                                GB32 = GB32 + (ZB^L)*(-1)^(L-1)/(L^1.5);
                            end
                            PSIA = 2*ET(i,m) - GA32*(R(i,m)/GA12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                            PSIB = 2*ET(i,m) - GB32*(R(i,m)/GB12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                            while (abs(ZA-ZB) > 0.00001)
                                
                                ZC = (ZA + ZB)/2;
                                GC12 = 0;
                                GC32 = 0;
                                GC52 = 0;
                                for L = 1:50
                                    GC12 = GC12 + (ZC^L)*(-1)^(L-1)/(L^0.5);
                                    GC32 = GC32 + (ZC^L)*(-1)^(L-1)/(L^1.5);
                                    GC52 = GC52 + (ZC^L)*(-1)^(L-1)/(L^2.5);
                                end
                                PSIC = 2*ET(i,m) - GC32*(R(i,m)/GC12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                                
                                if ((PSIA*PSIC) < 0)
                                    ZB = ZC;
                                    GB12 = GC12;
                                    GB32 = GC32;
                                    PSIB =PSIC;
                                else
                                    ZA = ZC;
                                    GA12 = GC12;
                                    GA32 = GC32;
                                    PSIA =PSIC;
                                end
                            end
                            Z(i,m) = ZC;
                            T(i,m) = R(i,m)^2 / (pi*GC12^2 );
                            P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                        end
                    end
                    
            end %if IT
            
        else
            F_new=F_new+ alpha(rk)*dt*F_tmp;
            %                 phi=phi+ alpha(rk)*dt*phi_t;
            %                 psi=psi+ alpha(rk)*dt*psi_t;
        end
    end  % RK
%%    
    F=F_new;
    
    for i=1:nx
        Mtemp=zeros(NV,pp);
        for K=1:NV
            Mtemp(K,:)=F(K,i,:);
        end
        F_loc(:,:)=Mtemp*Pleg;
        for m=1:pp
            SR(i,:) = wp * F_loc;
            SU(i,m) = sum(wp.*F_loc(:,m)'.* V);
            SE(i,m) = sum(wp.*F_loc(:,m)'.* V.^2)/2;
            SAV(i,m)= sum(wp.*F_loc(:,m)'.* abs(V));
            %         for K = 1: NV
            %             SR(i,m) = SR(i,m) + wp(K) * Floc(K,i,m);
            %             SU(i,m) = SU(i,m) + wp(K) *  Floc(K,i,m) * V(K);
            %             SE(i,m) = SE(i,m) + wp(K) *  Floc(K,i,m) * (0.5 * V(K) * V(K));
            %             SAV(i,m) = SAV(i,m) + wp(K) *  Floc(K,i,m) * abs(V(K));
            %         end
            R(i,m)    = SR(i,m);
            U(i,m)    = SU(i,m)/SR(i,m);
            ET(i,m)   = SE(i,m);
            AV(i,m)   = SAV(i,m);
        end
    end
    switch IT
        case 0
            %     if (IT == 0)
            for i=1:nx
                for m=1:pp
                    T(i,m)    = 4*ET(i,m)/R(i,m) - 2*U(i,m)^2;
                    Z(i,m)    = R(i,m) / SQRT(pi* T(i,m));
                    P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                end
            end
        case -1 %(IT == -1)
            for i=1:nx
                for m=1:pp
                    ZA = 0.0001;
                    ZB = 0.99;
                    GA12 = 0;
                    GB12 = 0;
                    GA32 = 0;
                    GB32 = 0;
                    for L = 1:50
                        GA12 = GA12 + (ZA^L)/(L^0.5);
                        GB12 = GB12 + (ZB^L)/(L^0.5);
                        GA32 = GA32 + (ZA^L)/(L^1.5);
                        GB32 = GB32 + (ZB^L)/(L^1.5);
                    end
                    PSIA = 2*ET(i,m) - GA32*(R(i,m)/GA12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                    PSIB = 2*ET(i,m) - GB32*(R(i,m)/GB12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                    while (abs(ZA-ZB) > 0.00001)
                        
                        ZC = (ZA + ZB)/2;
                        GC12 = 0;
                        GC32 = 0;
                        GC52 = 0;
                        for L = 1:50
                            GC12 = GC12 + (ZC^L)/(L^0.5);
                            GC32 = GC32 + (ZC^L)/(L^1.5);
                            GC52 = GC52 + (ZC^L)/(L^2.5);
                        end
                        PSIC = 2*ET(i,m) - GC32*(R(i,m)/GC12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                        
                        if ((PSIA*PSIC) < 0)
                            ZB = ZC;
                            GB12 = GC12;
                            GB32 = GC32;
                            PSIB =PSIC;
                        else
                            ZA = ZC;
                            GA12 = GC12;
                            GA32 = GC32;
                            PSIA =PSIC;
                        end
                    end
                    Z(i,m) = ZC;
                    T(i,m) = R(i,m)^2 / (pi*GC12^2 );
                    P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                end
            end
        case 1            %(IT == 1)
            for i=1:nx
                for m=1:pp
                    ZA = 0.0001;
                    ZB = 0.99;
                    GA12 = 0;
                    GB12 = 0;
                    GA32 = 0;
                    GB32 = 0;
                    for L = 1:50
                        GA12 = GA12 + (ZA^L)*(-1)^(L-1)/(L^0.5);
                        GB12 = GB12 + (ZB^L)*(-1)^(L-1)/(L^0.5);
                        GA32 = GA32 + (ZA^L)*(-1)^(L-1)/(L^1.5);
                        GB32 = GB32 + (ZB^L)*(-1)^(L-1)/(L^1.5);
                    end
                    PSIA = 2*ET(i,m) - GA32*(R(i,m)/GA12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                    PSIB = 2*ET(i,m) - GB32*(R(i,m)/GB12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                    while (abs(ZA-ZB) > 0.00001)
                        
                        ZC = (ZA + ZB)/2;
                        GC12 = 0;
                        GC32 = 0;
                        GC52 = 0;
                        for L = 1:50
                            GC12 = GC12 + (ZC^L)*(-1)^(L-1)/(L^0.5);
                            GC32 = GC32 + (ZC^L)*(-1)^(L-1)/(L^1.5);
                            GC52 = GC52 + (ZC^L)*(-1)^(L-1)/(L^2.5);
                        end
                        PSIC = 2*ET(i,m) - GC32*(R(i,m)/GC12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                        
                        if ((PSIA*PSIC) < 0)
                            ZB = ZC;
                            GB12 = GC12;
                            GB32 = GC32;
                            PSIB =PSIC;
                        else
                            ZA = ZC;
                            GA12 = GC12;
                            GA32 = GC32;
                            PSIA =PSIC;
                        end
                    end
                    Z(i,m) = ZC;
                    T(i,m) = R(i,m)^2 / (pi*GC12^2 );
                    P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
                end
            end
            
    end %if IT
 %%   
    %     for i=1:nx
    %         fp_plot((i-1)*pp+1:i*pp)=Floc(Kplot,i,1:pp);
    %         fm_plot((i-1)*pp+1:i*pp)=Floc(NV-Kplot+1,i,1:pp);
    %     end
    r_plot=reshape(R',nx*pp,1);
    u_plot=reshape(U',nx*pp,1);
    if bb==1
        title(sprintf('Time %8.4f',TIME));
        set(wave_handleu,'YData',r_plot); drawnow
        %         set(wave_handlev,'YData',u_plot); drawnow
    end
    %     if bb==1
    %         figure(1)
    %         title(sprintf('Iteration %d',ITER));
    %         set(wave_handleu,'YData',fp_plot); drawnow
    %         set(wave_handlev,'YData',fm_plot); drawnow
    %     end
    
    %fprintf('1X ELAPSED TIME: %f7.4,4 DENSITY AT X=4.0,Y=5.: %f7.4\n', TIME, R(NXP1/2))
    
    ITER = ITER + 1;
    
end
%% Enb for the program
% stop measuring time
toc

% 2000 CONTINUE
%     DO J = 1, NXP1
%     WRITE (10,*) X(J), R(J), P(J), Z(J), T(J), VIS(J)
%     END DO
% STOP
% END PROGRAM