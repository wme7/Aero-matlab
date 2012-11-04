% Arrays for Computation (Speed Up)
A=zeros(pp,pp);
b=zeros(1,pp);
c=zeros(1,pp);


MF=zeros(NV,NV);
F=zeros(NV,nx,pp);
FEQ=zeros(NV,nx,pp);
%Fd=zeros(NV,nx,pp);
F_tmp=zeros(NV,nx,pp);
F_new=zeros(NV,nx,pp);
FS=zeros(NV,nx,pp);
% Floc=zeros(NV,nx,pp);
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

FR=zeros(pp,1);
FU=zeros(pp,1);
FC=zeros(pp,1);
FN=zeros(pp,1);

% For Implicit-Explicit RK
F_s=zeros(NV,nx,pp,stage);
F_ns=zeros(NV,nx,pp,stage);

% For Plot Purpose
No=max(pp+1,12);
Fo=zeros(NV,No);
SRo=zeros(nx,No);
SUo=zeros(nx,No);
SEo=zeros(nx,No);
SAVo=zeros(nx,No);
Ro=zeros(nx,No);
Po=zeros(nx,No);
Uo=zeros(nx,No);
To=zeros(nx,No);
Zo=zeros(nx,No);
ETo=zeros(nx,No);
AVo=zeros(nx,No);
xo=zeros(nx,No);