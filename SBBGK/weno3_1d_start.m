%% SBBGK using WENO3
clc; clear all; close all;

%% Start
NX          = 101;
NXP1        = NX;
GHNC        = 1;
CFL         = 0.8;
OUTTIME     = 0.1;
IT       = 0;
if GHNC == 0 %THEN
    NV = 20;
else
    NV = 201;
end
%PI = ATAN2(1.;1.)*4.
GH= [-5.38748089001 -4.60368244955 -3.94476404012 -3.34785456738 ...
    -2.78880605843 -2.25497400209 -1.73853771212 -1.2340762154 ...
    -0.737473728545 -0.245340708301 0.245340708301 0.737473728545 ...
    1.2340762154 1.73853771212 2.25497400209 2.78880605843 ...
    3.34785456738 3.94476404012 4.60368244955 5.38748089001];
W = [0.898591961453 0.704332961176 0.62227869619 0.575262442852 ...
    0.544851742366 0.524080350949 0.509679027117 0.499920871336 ...
    0.493843385272 0.490921500667 0.490921500667 0.493843385272 ...
    0.499920871336 0.509679027117 0.524080350949 0.544851742366 ...
    0.575262442852 0.62227869619 0.704332961176 0.898591961453];
if GHNC == 0 %GO TO 10
    C = W;
    V = GH;
else
    V1  = -20;
    V2  = 20;
    DV  = (V2-V1)/(NV-1);
    for I = 1:NV
        V(I) = V1 + DV*(I-1);
    end
    for I = 2:NV-1
        C(I) = 64./45. * DV;
        
        if (mod(I,4)== 1)
            C(I) = 28./45. * DV;
        end
        
        if (mod(I,4)== 3)
            C(I) = 24./45. * DV;
        end
        
        C(1)    = 14./45. * DV;
        C(NV)   = C(1);
    end
end
        
%% Load IC
UL  = 0.;
TL  = 4.38385;
ZL  = 0.2253353;

UR  = 0.;
TR  = 8.972544;
ZR  = 0.1204582;
        
% create X domain        
DX = 1/(NX-1);
X(1)  = -0.5*DX;
for J = 2:NXP1
    X(J) = X(J-1) + DX;
end

%% load conditions into domain
for J = 1:NXP1
    if X(J)<= 0.5
        U(J) = UL;
        T(J) = TL;
        Z(J) = ZL;
    else
        U(J) = UR;
        T(J) = TR;
        Z(J) = ZR;
    end
    %evaluate equilibrium
    for K = 1:NV
        PP       = (V(K)-U(J))*(V(K)-U(J))/T(J);
        f(K,J)   = 1/((exp(PP)/Z(J)) + IT);
    end
end

%% Macrosopic Moments
for J = 1:NXP1
    SR = 0;
    SU = 0;
    SE = 0;
    SAV= 0;
    for K = 1:NV
        SR = SR + C(K) * f(K,J);
        SU = SU + C(K) * f(K,J) * V(K);
        SE = SE + C(K) * f(K,J) * (0.5 * V(K) * V(K));
        SAV = SAV + C(K) * f(K,J) * abs(V(K));
    end
    R(J)    = SR;
    j_x(J)  = SU;
    U(J)    = SU/SR;
    ET(J)   = SE;
    AV(J)   = SAV;
end
ITER  = 1;
TIME  = 0;
ISTOP = 0;

%% 1000 CONTINUE
for J = 1:NXP1
    c_kn = 0.001;
    tau = 0.001;
    
    vis(J) = c_kn / av(J);
end

DT = DX * CFL/V(NV);
TIME = TIME + DT;
DTDX = DT/DX;

if (TIME > OUTTIME) %THEN
    DTCFL = OUTTIME - (TIME - DTCFL);
    TIME = OUTTIME;
    DT = DTCFL/CFL;
    DTDX = DTCFL / DX;
    ISTOP = 1;
end
for J = 1:NXP1
    for K = 1:NV
        PP      = (V(K)-U(J)) * (V(K)-U(J)) /T(J);
        %!Feq(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
        F(K,J)   = 1/((EXP(PP)/Z(J)) + IT );
    end
end