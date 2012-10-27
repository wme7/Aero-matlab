%% Fugacity finder 1D
clc; clear all;

%% INPUTS
% starting point
et = 2.4;
r  = 1.;
u  = 2.;

% Statistic
IT = 1;
% tolerance
tol = 1e-7;

%% Bisection Method
za = 0.001;
zb = 0.999;
while abs(za-zb) >= tol
    ga12 = 0;
    gb12 = 0;
    ga32 = 0;
    gb32 = 0;
    for l = 1:50
        if IT == 1 %then
            ga12 = ga12 + (za^l)*(-1)^(l-1)/(l^0.5);
            gb12 = gb12 + (zb^l)*(-1)^(l-1)/(l^0.5);
            ga32 = ga32 + (za^l)*(-1)^(l-1)/(l^1.5);
            gb32 = gb32 + (zb^l)*(-1)^(l-1)/(l^1.5);
        else
            ga12 = ga12 + (za^l)/(l^0.5);
            gb12 = gb12 + (zb^l)/(l^0.5);
            ga32 = ga32 + (za^l)/(l^1.5);
            gb32 = gb32 + (zb^l)/(l^1.5);
        end
    end
    psia = 2*et - ga32*(r/ga12)^3/(2*pi) - r*u^2;
    psib = 2*et - gb32*(r/gb12)^3/(2*pi) - r*u^2;
    zc = (za + zb)/2;
    gc12 = 0;
    gc32 = 0;
    %gc52 = 0;
    for l = 1:50
        if IT==1 %then
            gc12 = gc12 + (zc^l)*(-1)^(l-1)/(l^0.5);
            gc32 = gc32 + (zc^l)*(-1)^(l-1)/(l^1.5);
            %gc52 = gc52 + (zc^l)*(-1)^(l-1)/(l^2.5);
        else
            gc12 = gc12 + (zc^l)/(l^0.5);
            gc32 = gc32 + (zc^l)/(l^1.5);
            %gc52 = gc52 + (zc^l)/(l^2.5);
        end
    end
    
    psic = 2*et - gc32*(r/gc12)^3/(2*pi) - r*u^2;
    
    if ((psia*psic) <= 0) %then
        zb = zc;
    else
        za = zc;
    end
end
%% OUTPUT
Z1 = zc;


%% Apply bisection method to the approx FD distribution Eq.
% Using same inputs
n = r;
E = et;
% u = u;

% FD case
r_a = 0.001; r_b = 0.999; tol = 1e-7;

for i = 1:1
    psi = @(r_x) 2*E(i)- FD(r_x,1.5)*(n(i)/FD(r_x,0.5))^3/(2*pi) ...
        - n(i)*(u(i)^2);
    r_p = bisection(psi,r_a,r_b,tol);
    x(i) = r_p;
    %t(i) = n(i)/(pi*FD(r_p,1));
    %p(i) = E(i) - 1/2*n(i)*(u(i)^2);
end

fprintf('Z1 = %2.12f\n',Z1);
fprintf('x  = %2.12f\n',x);