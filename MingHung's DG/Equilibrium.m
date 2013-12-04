for i=1:nx
    xi=(2*i-1)*dx/2; %evaluating the function `func' at the quadrature points
    x((i-1)*pp+1:i*pp)=xi+xl*dx/2;
    x_p(i,:)=xi+xl*dx/2;
    if(xi+xl(2)*dx/2 <= 0.5)
        U(i,:) = UL;
        T(i,:) = TL;
        Z(i,:) = ZL;
    else
        U(i,:) = UR;
        T(i,:) = TR;
        Z(i,:) = ZR;
    end
    
    for K = 1: NV
        for m=1:pp   %evaluating the function `func' at the quadrature points
            ffunc(m)  = 1/((exp((V(K)-U(i,m))^2/T(i,m))/Z(i,m)) + IT);
        end
        for j=0:p
            F(K,i,j+1)= sum (ffunc.*Pleg(j+1,:).*w)*(2*j+1)/2;
        end
    end
end
%Tmin=min(TR,TL);
%Tmax=max(TR,TL);

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
        
        R(i,m)    = SR(i,m);
        U(i,m)    = SU(i,m)/SR(i,m);
        ET(i,m)   = SE(i,m);
        AV(i,m)   = SAV(i,m);
        T(i,m)    = 4*ET(i,m)/R(i,m) - 2*U(i,m)^2;
        Z(i,m)    = R(i,m) / sqrt(pi* T(i,m));
    end
end