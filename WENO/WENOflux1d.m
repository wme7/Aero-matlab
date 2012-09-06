function [F_l,F_r] = WENOflux1d(u,a)
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
nx = length(u);
x = 4:nx-3;

% Flux Spliting Parameters
ap = max(a,0);
am = min(a,0);

% Initialize Arrays
flm = zeros(1,nx); flp = zeros(1,nx); frm = zeros(1,nx); frp = zeros(1,nx);

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
for i = x
    sl0mx = am*(13./12.)*(u(i-2)-2*u(i-1)+u(i))^2 ...
        + am*(1./4.)*(u(i-2)-4.*u(i-1)+3.*u(i))^2;
    sl1mx = am*(13./12.)*(u(i-1)-2*u(i)+u(i+1))^2 ...
        + am*(1./4.)*(u(i-1)-u(i+1))^2;
    sl2mx = am*(13./12.)*(u(i)-2*u(i+1)+u(i+2))^2 ...
        + am*(1./4.)*(3.*u(i)-4.*u(i+1)+u(i+2))^2;
    sl0px = ap*(13./12.)*(u(i-3)-2*u(i-2)+u(i-1))^2 ...
        + ap*(1./4.)*(u(i-3)-4.*u(i-2)+3.*u(i-1))^2;
    sl1px = ap*(13./12.)*(u(i-2)-2*u(i-1)+u(i))^2 ...
        + ap*(1./4.)*(u(i-2)-u(i))^2;
    sl2px = ap*(13./12.)*(u(i-1)-2*u(i)+u(i+1))^2 ...
        + ap*(1./4.)*(3.*u(i-1)-4.*u(i)+u(i+1))^2;
    sr0mx = am*(13./12.)*(u(i-1)-2*u(i)+u(i+1))^2 ...
        + am*(1./4.)*(u(i-1)-4.*u(i)+3.*u(i+1))^2;
    sr1mx = am*(13./12.)*(u(i)-2*u(i+1)+u(i+2))^2 ...
        + am*(1./4.)*(u(i)-u(i+2))^2;
    sr2mx = am*(13./12.)*(u(i+1)-2*u(i+2)+u(i+3))^2 ...
        + am*(1./4.)*(3.*u(i+1)-4.*u(i+2)+u(i+3))^2;
    sr0px = ap*(13./12.)*(u(i-2)-2*u(i-1)+u(i))^2 ...
        + ap*(1./4.)*(u(i-2)-4.*u(i-1)+3.*u(i))^2;
    sr1px = ap*(13./12.)*(u(i-1)-2*u(i)+u(i+1))^2 ...
        + ap*(1./4.)*(u(i-1)-u(i+1))^2;
    sr2px = ap*(13./12.)*(u(i)-2*u(i+1)+u(i+2))^2 ...
        + ap*(1./4.)*(3.*u(i)-4.*u(i+1)+u(i+2))^2;
    
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
    
    % negative & positive fluxes x
    flm(i) = am*(wl0mx*(c(1)*u(i-2)+c(3)*u(i-1)+c(2)*u(i))...
        +wl1mx*(c(2)*u(i-1)+c(3)*u(i)+c(1)*u(i+1))...
        +wl2mx*(c(5)*u(i)+c(4)*u(i+1)+c(2)*u(i+2)));
    flp(i) = ap*(wl0px*(c(2)*u(i-3)+c(4)*u(i-2)+c(5)*u(i-1))...
        +wl1px*(c(1)*u(i-2)+c(3)*u(i-1)+c(2)*u(i))...
        +wl2px*(c(2)*u(i-1)+c(3)*u(i)+c(1)*u(i+1)));
    frm(i) = am*(wr0mx*(c(1)*u(i-1)+c(3)*u(i)+c(2)*u(i+1))...
        +wr1mx*(c(2)*u(i)+c(3)*u(i+1)+c(1)*u(i+2))...
        +wr2mx*(c(5)*u(i+1)+c(4)*u(i+2)+c(2)*u(i+3)));
    frp(i) = ap*(wr0px*(c(2)*u(i-2)+c(4)*u(i-1)+c(5)*u(i))...
        +wr1px*(c(1)*u(i-1)+c(3)*u(i)+c(2)*u(i+1))...
        +wr2px*(c(2)*u(i)+c(3)*u(i+1)+c(1)*u(i+2)));
end

% fluxes x
F_l = flm + flp;
F_r = frm + frp;