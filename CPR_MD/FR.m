% Using Flux Reconstruction Procedure

%% Simulation Parameters
a = -1; % advection speed
%cfl = 0.2; % CFL condition
%tEnd = 1; % final time
np = 4; % degree of accuaracy
nc = 10; % number of elements
nf = nc+1;

%% Prepocessing 
% Define our Flux function
flux = @(w) a*w; %w.^2/2; 

% Derivate of the flux function
dflux = @(w) a*ones(size(w)); %w;

% Build 1d mesh
grid = mesh1d([0 1],nc,'LGL',np-1);
dx = grid.elementSize; Ja = grid.Jacobian; x = grid.nodeCoordinates;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',np); % g: one-order higher
dRR = RR.eval_dP(grid.solutionPoints); dRL = -flipud(dRR);

% Build Lagrange k-Polynomials
L = LagrangePolynomial(grid.solutionPoints);

% IC
u0 = IC(x,1);

QC = u0;

Coe_Lag = zeros(np,2);
Coe_Lag(:,1) = double(subs(L.lagrangePolynomial,-1));
Coe_Lag(:,2) = double(subs(L.lagrangePolynomial,1));
Coe_DLag = double(subs(L.dlagrangePolynomial,grid.solutionPoints));

%************* Reconstruct all fluxes from interpolation******************!
for ic=1:nc
    for i=1:np	%every interface(FPs)
        Flux(i,ic) = a * QC(i,ic) / Ja;
        if(i==1) %then	%store the SV interface valuces
            QcSVB(1,ic) = 0;
            for j=1:np
                QcSVB(1,ic) = QcSVB(1,ic) + Coe_Lag(j,1) * QC(j,ic);
            end
        elseif (i==np) %then
            QcSVB(2,ic+1) = 0;
            for j=1:np
                QcSVB(2,ic+1) = QcSVB(2,ic+1) + Coe_Lag(j,2) * QC(j,ic);
            end
        end
    end
    Flux_bond(1,ic) = 0;
    Flux_bond(2,ic) = 0;
    for j=1:np
        Flux_bond(1,ic) = Flux_bond(1,ic) + Coe_Lag(j,1) * Flux(j,ic);
        Flux_bond(2,ic) = Flux_bond(2,ic) + Coe_Lag(j,2) * Flux(j,ic);
    end
end

%Flux_bond
%f_bd = flux(u0);
%[f_bd(1,:);f_bd(np,:)]/Ja

    
%************Use Riemman solver to reconstruct SV boundary fluxes*********!
for i=1:nf
    if(i==1) %then	%left boundary, periodic
        QcSVB(2,i)=0;
        for j=1:np
            QcSVB(2,i) = QcSVB(2,i) + Coe_Lag(j,2) * QC(j,nc);
        end
    elseif(i==nf) %then	%right boundary
        QcSVB(1,i)=0;
        for j=1:np
            QcSVB(1,i) = QcSVB(1,i) + Coe_Lag(j,1) * QC(j,1);
        end
    end
    
    QCL=QcSVB(2,i);
    QCR=QcSVB(1,i);
    
    
    Cs=abs(a);
    FluxSVB(i) = 0.5d0*(QCL*a+QCR*a-Cs*(QCR-QCL)); %ok
end
    
%Flux_com(1:nf) = FluxSVB(1:nf) / Ja;
Flux_com = FluxSVB / Ja;

for ic=1:nc
    Flux_cor_l(ic)=Flux_com(ic)-Flux_bond(1,ic);
    Flux_cor_r(ic)=Flux_com(ic+1)-Flux_bond(2,ic);
end

for ic=1:nc
    for i=1:np
        DFlux(i,ic)=0;
        for j=1:np
            DFlux(i,ic)=DFlux(i,ic)+Coe_DLag(j,i)*Flux(j,ic);
        end
    end
end
DFlux;
% 
% for ic=1,nc	%!update CVs using correction via flux reconstructions
%     for i=1,np
%         Res(i,ic)=DFlux(i,ic) + Flux_cor_l(ic)*Coe_RadauR(i) + Flux_cor_r(ic)*Coe_RadauL(i);
%     end
% end
% 
% Res = -Res;