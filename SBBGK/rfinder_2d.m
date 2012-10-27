%% Fugacity finder 2D
clc; clear all;

%% INPUTS
% starting point
ET = 0.22263;
R  = 0.521295;
UX = -0.42308;
UY = 0.;

% Statistic
IT = 1;
% tolerance
tol = 1e-7;

%% Bisection Method
ZA = 0.001;
ZB = 0.999;
while abs(ZA-ZB) >= tol
    GA1 = 0;
    GB1 = 0;
    GA2 = 0;
    GB2 = 0;
    for L = 1:50
        GA1 = GA1 + (ZA^L) * (-1)^(L-1)/L;
        GB1 = GB1 + (ZB^L) * (-1)^(L-1)/L;
        GA2 = GA2 + (ZA^L) * (-1)^(L-1)/(L^2);
        GB2 = GB2 + (ZB^L) * (-1)^(L-1)/(L^2);
    end
    PSIA = 2*ET - (GA2*(R/GA1)^2)/pi - R*(UX*UX+UY*UY);
    PSIB = 2*ET - (GB2*(R/GB1)^2)/pi - R*(UX*UX+UY*UY);
    ZC = (ZA+ZB)/2;
    GC1 = 0;
    GC2 = 0;
    for L = 1:50
        GC1 = GC1 + (ZC^L) * (-1)^(L-1)/L;
        GC2 = GC2 + (ZC^L) * (-1)^(L-1)/(L^2);
    end
    PSIC = 2*ET - (GC2*(R/GC1)^2)/pi - R*(UX*UX+UY*UY);
    if PSIA*PSIC <= 0 %THEN
        ZB = ZC;
    else
        ZA = ZC;
    end
end

%% OUTPUT
Z1 = ZC;

%% Apply bisection method to the approx FD distribution Eq.

% Using same inputs
n = R;
u1= UX;
u2= UY;
E = ET;

% FD case
r_a = 0.001; r_b = 0.999;

for j = 1:1
    for i = 1:1
        psi = @(r_x)2*E(j,i) - FD(r_x,2)/pi*(n(j,i)/FD(r_x,1)).^2 ...
            - n(j,i).*(u1(j,i).^2 + u2(j,i).^2);
        r_p = bisection(psi,r_a,r_b,tol);
        r(j,i) = r_p;
        %t(j,i) = n(i)/(pi*FermiDiract(r_p,1));
        %p(j,i) = E(j,i) - 1/2*n(j,i)*(u1(j,i).^2 + u2(j,i).^2);
    end
end

fprintf('Z1 = %2.12f\n',Z1);
fprintf('r  = %2.12f\n',r);

% if Z1 == r
%     fprintf('congratz!\n')
% else
%     error('failed')
% end