function ft=BGKimL(u)
% Where f is the right hand side of the steady-state example
%
%
global A b c Pleg w wp NV nx p dx dt IT BC_type V VIS F gamma

pp=p+1;
NVh=NV/2;

Fv=reshape(u,NV,nx,pp);

fta=zeros(NV,nx,pp);
FS=zeros(NV,nx,pp);
FEQ=zeros(NV,nx,pp);
FC=zeros(pp,1);
FB=zeros(pp,1);
FR=zeros(pp,1);
F_s=Fv;

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

for i=1:nx
    Mtemp=zeros(NV,pp);
    for K=1:NV
        Mtemp(K,:)=Fv(K,i,:);
    end
    F_loc(:,:)=Mtemp*Pleg;
    for m=1:pp
        SR(i,:) = wp * F_loc;
        SU(i,m) = sum(wp.*F_loc(:,m)'.* V);
        SE(i,m) = sum(wp.*F_loc(:,m)'.* V.^2)/2;
        SAV(i,m)= sum(wp.*F_loc(:,m)'.* abs(V));
        
        R(i,m)    = SR(i,m);
        U(i,m)    = SU(i,m)/SR(i,m);
        ET(i,m)   = SE(i,m);
        AV(i,m)   = SAV(i,m);
    end
end

if (IT == 0)
    for i=1:nx
        for m=1:pp
            T(i,m)    = 4*ET(i,m)/R(i,m) - 2*U(i,m)^2;
            Z(i,m)    = R(i,m) / sqrt(pi* T(i,m));
            P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
        end
    end
else
    for i=1:nx
        for m=1:pp
            ZA = 0.0001;
            ZB = 0.99;
            while (abs(ZA-ZB) > 0.00001)
                GA12 = 0;
                GB12 = 0;
                GA32 = 0;
                GB32 = 0;
                for L = 1:50
                    if (IT == 1)
                        GA12 = GA12 + (ZA^L)*(-1)^(L-1)/(L^0.5);
                        GB12 = GB12 + (ZB^L)*(-1)^(L-1)/(L^0.5);
                        GA32 = GA32 + (ZA^L)*(-1)^(L-1)/(L^1.5);
                        GB32 = GB32 + (ZB^L)*(-1)^(L-1)/(L^1.5);
                    else
                        GA12 = GA12 + (ZA^L)/(L^0.5);
                        GB12 = GB12 + (ZB^L)/(L^0.5);
                        GA32 = GA32 + (ZA^L)/(L^1.5);
                        GB32 = GB32 + (ZB^L)/(L^1.5);
                    end
                end
                PSIA = 2*ET(i,m) - GA32*(R(i,m)/GA12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                PSIB = 2*ET(i,m) - GB32*(R(i,m)/GB12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                ZC = (ZA + ZB)/2;
                GC12 = 0;
                GC32 = 0;
                GC52 = 0;
                for L = 1:50
                    if  (IT == 1)
                        GC12 = GC12 + (ZC^L)*(-1)^(L-1)/(L^0.5);
                        GC32 = GC32 + (ZC^L)*(-1)^(L-1)/(L^1.5);
                        GC52 = GC52 + (ZC^L)*(-1)^(L-1)/(L^2.5);
                    else
                        GC12 = GC12 + (ZC^L)/(L^0.5);
                        GC32 = GC32 + (ZC^L)/(L^1.5);
                        GC52 = GC52 + (ZC^L)/(L^2.5);
                    end
                end
                PSIC = 2*ET(i,m) - GC32*(R(i,m)/GC12)^3/(2*pi) - R(i,m)*U(i,m)^2;
                
                if ((PSIA*PSIC) < 0)
                    ZB = ZC;
                else
                    ZA = ZC;
                end
            end
            Z(i,m) = ZC;
            T(i,m) = R(i,m)^2 / (pi*GC12^2 );
            P(i,m) = ET(i,m) - 0.5 * R(i,m) * U(i,m)^2;
            
        end
    end
end %if IT



for i = 1:nx
    for K = 1:NV
        for m=1:pp
            FEQ(K,i,m)   = 1/((exp( (V(K)-U(i,m))^2 /T(i,m))/Z(i,m)) + IT );
        end
    end
end

for i=1:nx
    for K = 1: NV
        FC(:)=Fv(K,i,:);
        FB(:)=FEQ(K,i,:);
        FC=(FC'*Pleg-FB')';        
        for j=0:p
            %         FS(K,i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*(2*j+1)/2/VIS(i,j+1);
            FS(K,i,j+1)=sum (FC'.*Pleg(j+1,:).*w)*dx/2/VIS(i,j+1);
        end
    end
    
end
for K=1:NVh
    if BC_type == 0
        %BC no-flux
        FC(:)=Fv(K,1,:);
        FU=FC;
        FR(:)=FS(K,1,:);
        F_s(K,1,:)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c'-FR)' .* b);
    elseif BC_type == -1            %BC reflecting
        FC(:)=Fv(K,1,:);
        FU(:)=Fv(NV-K+1,1,:);
        FR(:)=FS(K,1,:);
        F_s(K,1,:)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c'-FR)' .* b);
    else
    end
    for i=2:nx
        FU(:)=Fv(K,i-1,:);
        FC(:)=Fv(K,i,:);
        FR(:)=FS(K,i,:);
        F_s(K,i,:)=( (V(K)*A'*FC -V(K)* sum(FC) +V(K)* sum(FU) * c'-FR)' .* b);
    end
    
    for i=1:nx-1
        FU(:)=Fv(NVh+K,i+1,:);
        FC(:)=Fv(NVh+K,i,:);
        FR(:)=FS(NVh+K,i,:);
        F_s(NVh+K,i,:)=( (-V(NVh+K)*A'*(-FC) -V(NVh+K)* sum(FU'.* c) + V(NVh+K)*sum(FC'.* c) * c'-FR)' .* b);
    end
    if BC_type == 0
        %BC no-flux
        FC(:)=Fv(NVh+K,nx,:);
        FU=FC;
        FR(:)=FS(NVh+K,nx,:);
        F_s(NVh+K,nx,:)=( (-V(NVh+K)*A'*(-FC) -V(NVh+K)* sum(FU) +V(NVh+K)* sum(FC'.*c) * c'-FR)' .* b);
    elseif BC_type == -1
        %BC reflexting
        FC(:)=Fv(NVh+K,nx,:);
        FU(:)=Fv(NVh-K+1,nx,:);
        FR(:)=FS(NVh+K,nx,:);
        F_s(NVh+K,nx,:)=( (-V(NVh+K)*A'*(-FC) -V(NVh+K)* sum(FU) +V(NVh+K)* sum(FC'.*c) * c'-FR)' .* b);
    end
end % loop for NV

fta= Fv - F - gamma*dt*F_s;
ft=reshape(fta,NV*nx*pp,1);