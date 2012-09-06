function [F_l,F_r,G_l,G_r] = WENOflux2d(u,a,b)
%% Weighted Essentially Non-Oscilatory for 5th order Accuracy
% WENO 3 O(h^5) FDM implementation subroutine for computing the numerical
% fluxes at the right and left boundaries of the every cell, that is:
%
% $v_{i+\frac{1}{2}}^{-}$ and % $v_{i-\frac{1}{2}}^{+}$. 
%
% This algorithm based on lectures notes of:
%
%   Chi-Wang Shu; High-Order ENO and WENO schemes for Computational
%   Fluid Dynamics, High-Order Methods for Computational Physics.
%   Springer 1999.
%
% The basic idea is the following: instead of using only one of ENO
% stencils to form a reconstruction, one uses a convex combination of all
% of them as:
%
% For the right boundary: 
% $v_{i+\frac{1}{2}}^{-} = \sum_{r=0}^{k-1}\omega_{r} v_{i+\frac{1}{2}}^{(r)}$
%
% For the left boundary: 
% $v_{i-\frac{1}{2}}^{+} = \sum_{r=0}^{k-1}\widetilde{\omega}_{r} v_{i-\frac{1}{2}}^{(r)}$

%% Parameters:
[ny,nx] = size(u);
x = 4:nx-3;
y = 4:ny-3;

% Flux Spliting Parameters
ap = max(a,0);
am = min(a,0);
bp = max(b,0);
bm = min(b,0);

% Initialize Arrays
flm = zeros(ny,nx); flp = zeros(ny,nx); frm = zeros(ny,nx); frp = zeros(ny,nx);
glm = zeros(ny,nx); glp = zeros(ny,nx); grm = zeros(ny,nx); grp = zeros(ny,nx);

% %% Polynomial Reconstruction Coeficients (C_rj)
% switch Degree
%     
%     case{1} % 1D Problem, % 1rd Order ENO, 3th order WENO
%         c = [ 1; ...
%               1 ];
%     case{2} % 2D Problem, % 2rd Order ENO, 4th order WENO
%         c = [ 3/2 -1/2; ...
%               1/2  1/2; ...
%              -1/2  3/2 ];
%     case{3} % 3D Problem, % 3rd Order ENO, 5th order WENO
%         c = [ 11/6 -7/6  1/3; ...
%                1/3  5/6 -1/6; ...
%               -1/6  5/6  1/3; ...
%                1/3 -7/6 11/6 ];
%     otherwise
%         error('only available cases: Degree = 1, 2 and 3')
% end

    c = [-1/6 1/3 5/6 -7/6 11/6];
    
%% Compute Weno Fluxes:
for j = y
    for i = x
        sl0my = bm*(13./12.)*(u(j-2,i)-2*u(j-1,i)+u(j,i))^2 ...
            + bm*(1./4.)*(u(j-2,i)-4.*u(j-1,i)+3.*u(j,i))^2;
        sl1my = bm*(13./12.)*(u(j-1,i)-2*u(j,i)+u(j+1,i))^2 ...
            + bm*(1./4.)*(u(j-1,i)-u(j+1,i))^2;
        sl2my = bm*(13./12.)*(u(j,i)-2*u(j+1,i)+u(j+2,i))^2 ...
            + bm*(1./4.)*(3.*u(j,i)-4.*u(j+1,i)+u(j+2,i))^2;
        sl0py = bp*(13./12.)*(u(j-3,i)-2*u(j-2,i)+u(j-1,i))^2 ...
            + bp*(1./4.)*(u(j-3,i)-4.*u(j-2,i)+3.*u(j-1,i))^2;
        sl1py = bp*(13./12.)*(u(j-2,i)-2*u(j-1,i)+u(j,i))^2 ...
            + bp*(1./4.)*(u(j-2,i)-u(j,i))^2;
        sl2py = bp*(13./12.)*(u(j-1,i)-2*u(j,i)+u(j+1,i))^2 ...
            + bp*(1./4.)*(3.*u(j-1,i)-4.*u(j,i)+u(j+1,i))^2;
        sr0my = bm*(13./12.)*(u(j-1,i)-2*u(j,i)+u(j+1,i))^2 ...
            + bm*(1./4.)*(u(j-1,i)-4.*u(j,i)+3.*u(j+1,i))^2;
        sr1my = bm*(13./12.)*(u(j,i)-2*u(j+1,i)+u(j+2,i))^2 ...
            + bm*(1./4.)*(u(j,i)-u(j+2,i))^2;
        sr2my = bm*(13./12.)*(u(j+1,i)-2*u(j+2,i)+u(j+3,i))^2 ...
            + bm*(1./4.)*(3.*u(j+1,i)-4.*u(j+2,i)+u(j+3,i))^2;
        sr0py = bp*(13./12.)*(u(j-2,i)-2*u(j-1,i)+u(j,i))^2 ...
            + bp*(1./4.)*(u(j-2,i)-4.*u(j-1,i)+3.*u(j,i))^2;
        sr1py = bp*(13./12.)*(u(j-1,i)-2*u(j,i)+u(j+1,i))^2 ...
            + bp*(1./4.)*(u(j-1,i)-u(j+1,i))^2;
        sr2py = bp*(13./12.)*(u(j,i)-2*u(j+1,i)+u(j+2,i))^2 ...
            + bp*(1./4.)*(3.*u(j,i)-4.*u(j+1,i)+u(j+2,i))^2;
        
        sl0mx = am*(13./12.)*(u(j,i-2)-2*u(j,i-1)+u(j,i))^2 ...
            + am*(1./4.)*(u(j,i-2)-4.*u(j,i-1)+3.*u(j,i))^2;
        sl1mx = am*(13./12.)*(u(j,i-1)-2*u(j,i)+u(j,i+1))^2 ...
            + am*(1./4.)*(u(j,i-1)-u(j,i+1))^2;
        sl2mx = am*(13./12.)*(u(j,i)-2*u(j,i+1)+u(j,i+2))^2 ...
            + am*(1./4.)*(3.*u(j,i)-4.*u(j,i+1)+u(j,i+2))^2;
        sl0px = ap*(13./12.)*(u(j,i-3)-2*u(j,i-2)+u(j,i-1))^2 ...
            + ap*(1./4.)*(u(j,i-3)-4.*u(j,i-2)+3.*u(j,i-1))^2;
        sl1px = ap*(13./12.)*(u(j,i-2)-2*u(j,i-1)+u(j,i))^2 ...
            + ap*(1./4.)*(u(j,i-2)-u(j,i))^2;
        sl2px = ap*(13./12.)*(u(j,i-1)-2*u(j,i)+u(j,i+1))^2 ...
            + ap*(1./4.)*(3.*u(j,i-1)-4.*u(j,i)+u(j,i+1))^2;
        sr0mx = am*(13./12.)*(u(j,i-1)-2*u(j,i)+u(j,i+1))^2 ...
            + am*(1./4.)*(u(j,i-1)-4.*u(j,i)+3.*u(j,i+1))^2;
        sr1mx = am*(13./12.)*(u(j,i)-2*u(j,i+1)+u(j,i+2))^2 ...
            + am*(1./4.)*(u(j,i)-u(j,i+2))^2;
        sr2mx = am*(13./12.)*(u(j,i+1)-2*u(j,i+2)+u(j,i+3))^2 ...
            + am*(1./4.)*(3.*u(j,i+1)-4.*u(j,i+2)+u(j,i+3))^2;
        sr0px = ap*(13./12.)*(u(j,i-2)-2*u(j,i-1)+u(j,i))^2 ...
            + ap*(1./4.)*(u(j,i-2)-4.*u(j,i-1)+3.*u(j,i))^2;
        sr1px = ap*(13./12.)*(u(j,i-1)-2*u(j,i)+u(j,i+1))^2 ...
            + ap*(1./4.)*(u(j,i-1)-u(j,i+1))^2;
        sr2px = ap*(13./12.)*(u(j,i)-2*u(j,i+1)+u(j,i+2))^2 ...
            + ap*(1./4.)*(3.*u(j,i)-4.*u(j,i+1)+u(j,i+2))^2;
        
        al0my  = 1. / (10. * (eps + sl0my))^2;
        al1my  = 6. / (10. * (eps + sl1my))^2;
        al2my  = 3. / (10. * (eps + sl2my))^2;
        al0py  = 1. / (10. * (eps + sl0py))^2;
        al1py  = 6. / (10. * (eps + sl1py))^2;
        al2py  = 3. / (10. * (eps + sl2py))^2;
        ar0my  = 3. / (10. * (eps + sr0my))^2;
        ar1my  = 6. / (10. * (eps + sr1my))^2;
        ar2my  = 1. / (10. * (eps + sr2my))^2;
        ar0py  = 3. / (10. * (eps + sr0py))^2;
        ar1py  = 6. / (10. * (eps + sr1py))^2;
        ar2py  = 1. / (10. * (eps + sr2py))^2;
        
        al0mx  = 1. / (10. * (eps + sl0mx))^2;
        al1mx  = 6. / (10. * (eps + sl1mx))^2;
        al2mx  = 3. / (10. * (eps + sl2mx))^2;
        al0px  = 1. / (10. * (eps + sl0px))^2;
        al1px  = 6. / (10. * (eps + sl1px))^2;
        al2px  = 3. / (10. * (eps + sl2px))^2;
        ar0mx  = 3. / (10. * (eps + sr0mx))^2;
        ar1mx  = 6. / (10. * (eps + sr1mx))^2;
        ar2mx  = 1. / (10. * (eps + sr2mx))^2;
        ar0px  = 3. / (10. * (eps + sr0px))^2;
        ar1px  = 6. / (10. * (eps + sr1px))^2;
        ar2px  = 1. / (10. * (eps + sr2px))^2;
        
        % weightings y
        wl0my  = al0my / (al0my+al1my+al2my);
        wl1my  = al1my / (al0my+al1my+al2my);
        wl2my  = al2my / (al0my+al1my+al2my);
        wl0py  = al0py / (al0py+al1py+al2py);
        wl1py  = al1py / (al0py+al1py+al2py);
        wl2py  = al2py / (al0py+al1py+al2py);
        wr0my  = ar0my / (ar0my+ar1my+ar2my);
        wr1my  = ar1my / (ar0my+ar1my+ar2my);
        wr2my  = ar2my / (ar0my+ar1my+ar2my);
        wr0py  = ar0py / (ar0py+ar1py+ar2py);
        wr1py  = ar1py / (ar0py+ar1py+ar2py);
        wr2py  = ar2py / (ar0py+ar1py+ar2py);
        
        % weightings x
        wl0mx  = al0mx / (al0mx+al1mx+al2mx);
        wl1mx  = al1mx / (al0mx+al1mx+al2mx);
        wl2mx  = al2mx / (al0mx+al1mx+al2mx);
        wl0px  = al0px / (al0px+al1px+al2px);
        wl1px  = al1px / (al0px+al1px+al2px);
        wl2px  = al2px / (al0px+al1px+al2px);
        wr0mx  = ar0mx / (ar0mx+ar1mx+ar2mx);
        wr1mx  = ar1mx / (ar0mx+ar1mx+ar2mx);
        wr2mx  = ar2mx / (ar0mx+ar1mx+ar2mx);
        wr0px  = ar0px / (ar0px+ar1px+ar2px);
        wr1px  = ar1px / (ar0px+ar1px+ar2px);
        wr2px  = ar2px / (ar0px+ar1px+ar2px);
        
        % negative & positive fluxes y
        glm(j,i) = bm*(wl0my*(c(1)*u(j-2,i)+c(3)*u(j-1,i)+c(2)*u(j,i))...
            +wl1my*(c(2)*u(j-1,i)+c(3)*u(j,i)+c(1)*u(j+1,i))...
            +wl2my*(c(5)*u(j,i)+c(4)*u(j+1,i)+c(2)*u(j+2,i)));
        glp(j,i) = bp*(wl0py*(c(2)*u(j-3,i)+c(4)*u(j-2,i)+c(5)*u(j-1,i))...
            +wl1py*(c(1)*u(j-2,i)+c(3)*u(j-1,i)+c(2)*u(j,i))...
            +wl2py*(c(2)*u(j-1,i)+c(3)*u(j,i)+c(1)*u(j+1,i)));
        grm(j,i) = bm*(wr0my*(c(1)*u(j-1,i)+c(3)*u(j,i)+c(2)*u(j+1,i))...
            +wr1my*(c(2)*u(j,i)+c(3)*u(j+1,i)+c(1)*u(j+2,i))...
            +wr2my*(c(5)*u(j+1,i)+c(4)*u(j+2,i)+c(2)*u(j+3,i)));
        grp(j,i) = bp*(wr0py*(c(2)*u(j-2,i)+c(4)*u(j-1,i)+c(5)*u(j,i))...
            +wr1py*(c(1)*u(j-1,i)+c(3)*u(j,i)+c(2)*u(j+1,i))...
            +wr2py*(c(2)*u(j,i)+c(3)*u(j+1,i)+c(1)*u(j+2,i)));
        
        % negative & positive fluxes x
        flm(j,i) = am*(wl0mx*(c(1)*u(j,i-2)+c(3)*u(j,i-1)+c(2)*u(j,i))...
            +wl1mx*(c(2)*u(j,i-1)+c(3)*u(j,i)+c(1)*u(j,i+1))...
            +wl2mx*(c(5)*u(j,i)+c(4)*u(j,i+1)+c(2)*u(j,i+2)));
        flp(j,i) = ap*(wl0px*(c(2)*u(j,i-3)+c(4)*u(j,i-2)+c(5)*u(j,i-1))...
            +wl1px*(c(1)*u(j,i-2)+c(3)*u(j,i-1)+c(2)*u(j,i))...
            +wl2px*(c(2)*u(j,i-1)+c(3)*u(j,i)+c(1)*u(j,i+1)));
        frm(j,i) = am*(wr0mx*(c(1)*u(j,i-1)+c(3)*u(j,i)+c(2)*u(j,i+1))...
            +wr1mx*(c(2)*u(j,i)+c(3)*u(j,i+1)+c(1)*u(j,i+2))...
            +wr2mx*(c(5)*u(j,i+1)+c(4)*u(j,i+2)+c(2)*u(j,i+3)));
        frp(j,i) = ap*(wr0px*(c(2)*u(j,i-2)+c(4)*u(j,i-1)+c(5)*u(j,i))...
            +wr1px*(c(1)*u(j,i-1)+c(3)*u(j,i)+c(2)*u(j,i+1))...
            +wr2px*(c(2)*u(j,i)+c(3)*u(j,i+1)+c(1)*u(j,i+2)));    
    end
end

% fluxes x
F_l = flm + flp;
F_r = frm + frp;

% fluxes y
G_l = glm + glp;
G_r = grm + grp;