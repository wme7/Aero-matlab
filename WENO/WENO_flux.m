function [F_l,F_r,G_l,G_r] = WENO_flux(u,a,b,dx,dy,dt,w,Degree,d)
%% Weighted Essentially Non-Oscilatory for 5th order Accuracy
% WENO 5 implementation subroutine for computing numerical fluxes at the
% right and left boundaries of the every cell, that is 
% $v_{i+\frac{1}{2}}^{-}$ and % $v_{i-\frac{1}{2}}^{+}$. To solve for a 1D 
% scalar advection equation. 
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
% For the right boundary: $v_{i+\frac{1}{2}}^{-}=\sum_{r=0}^{k-1}\omega_{r}
% v_{i+\frac{1}{2}}^{(r)}$
%
% For the left boundary: $v_{i-\frac{1}{2}}^{+}=\sum_{r=0}^{k-1}
% \widetilde{\omega}_{r} v_{i-\frac{1}{2}}^{(r)}$
%% Compute WENO Fluxes:
[ny,nx] = size(u);
x = 4:nx-3;
y = 4:ny-3;

% Flux Spliting Parameters
vxp = max(a,0);
vxm = min(a,0);
vyp = max(b,0);
vym = min(b,0);

% Polynomial Reconstruction Coeficients (C_rj)
switch Degree
    
    case{1} % 1D Problem, % 1rd Order ENO, 3th order WENO
        c = [ 1;
              1 ];
    case{2} % 2D Problem, % 2rd Order ENO, 4th order WENO
        c = [ 3/2 -1/2;
              1/2  1/2;
             -1/2  3/2 ];
    case{3} % 3D Problem, % 3rd Order ENO, 5th order WENO
        c = [ 11/6 -7/6  1/3; ...
               1/3  5/6 -1/6; ...
              -1/6  5/6  1/3; ...
               1/3 -7/6 11/6 ];
    otherwise
        error('only available cases: Degree = 1, 2 and 3')
end
      
%% Reconstructions:

    % right    
        vr(i) = c_rj(1,1)*v() + c_rj(1,2)*v() + c_rj(1,3)*v()
    % left 
 

DO K = 1, NV
DO L = 1, NV
 
    
    DO I = 4, NX-3
    DO J = 4, NY-3     
             sl0mx = vxm*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 + vxm*(1./4.)*(f(k,l,i-2,j)-4.*f(k,l,i-1,j)+3.*f(k,l,i,j))**2
             sl1mx = vxm*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 + vxm*(1./4.)*(f(k,l,i-1,j)-f(k,l,i+1,j))**2
             sl2mx = vxm*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 + vxm*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i+1,j)+f(k,l,i+2,j))**2
             sl0px = vxp*(13./12.)*(f(k,l,i-3,j)-2*f(k,l,i-2,j)+f(k,l,i-1,j))**2 + vxp*(1./4.)*(f(k,l,i-3,j)-4.*f(k,l,i-2,j)+3.*f(k,l,i-1,j))**2
             sl1px = vxp*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 + vxp*(1./4.)*(f(k,l,i-2,j)-f(k,l,i,j))**2
             sl2px = vxp*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 + vxp*(1./4.)*(3.*f(k,l,i-1,j)-4.*f(k,l,i,j)+f(k,l,i+1,j))**2
             sr0mx = vxm*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 + vxm*(1./4.)*(f(k,l,i-1,j)-4.*f(k,l,i,j)+3.*f(k,l,i+1,j))**2
             sr1mx = vxm*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 + vxm*(1./4.)*(f(k,l,i,j)-f(k,l,i+2,j))**2
             sr2mx = vxm*(13./12.)*(f(k,l,i+1,j)-2*f(k,l,i+2,j)+f(k,l,i+3,j))**2 + vxm*(1./4.)*(3.*f(k,l,i+1,j)-4.*f(k,l,i+2,j)+f(k,l,i+3,j))**2
             sr0px = vxp*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 + vxp*(1./4.)*(f(k,l,i-2,j)-4.*f(k,l,i-1,j)+3.*f(k,l,i,j))**2
             sr1px = vxp*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 + vxp*(1./4.)*(f(k,l,i-1,j)-f(k,l,i+1,j))**2
             sr2px = vxp*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 + vxp*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i+1,j)+f(k,l,i+2,j))**2
             
             sl0my = vym*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 + vym*(1./4.)*(f(k,l,i,j-2)-4.*f(k,l,i,j-1)+3.*f(k,l,i,j))**2
             sl1my = vym*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 + vym*(1./4.)*(f(k,l,i,j-1)-f(k,l,i,j+1))**2
             sl2my = vym*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 + vym*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i,j+1)+f(k,l,i,j+2))**2
             sl0py = vyp*(13./12.)*(f(k,l,i,j-3)-2*f(k,l,i,j-2)+f(k,l,i,j-1))**2 + vyp*(1./4.)*(f(k,l,i,j-3)-4.*f(k,l,i,j-2)+3.*f(k,l,i,j-1))**2
             sl1py = vyp*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 + vyp*(1./4.)*(f(k,l,i,j-2)-f(k,l,i,j))**2
             sl2py = vyp*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 + vyp*(1./4.)*(3.*f(k,l,i,j-1)-4.*f(k,l,i,j)+f(k,l,i,j+1))**2
             sr0my = vym*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 + vym*(1./4.)*(f(k,l,i,j-1)-4.*f(k,l,i,j)+3.*f(k,l,i,j+1))**2
             sr1my = vym*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 + vym*(1./4.)*(f(k,l,i,j)-f(k,l,i,j+2))**2
             sr2my = vym*(13./12.)*(f(k,l,i,j+1)-2*f(k,l,i,j+2)+f(k,l,i,j+3))**2 + vym*(1./4.)*(3.*f(k,l,i,j+1)-4.*f(k,l,i,j+2)+f(k,l,i,j+3))**2
             sr0py = vyp*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 + vyp*(1./4.)*(f(k,l,i,j-2)-4.*f(k,l,i,j-1)+3.*f(k,l,i,j))**2
             sr1py = vyp*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 + vyp*(1./4.)*(f(k,l,i,j-1)-f(k,l,i,j+1))**2
             sr2py = vyp*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 + vyp*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i,j+1)+f(k,l,i,j+2))**2
             
             al0mx  = 1. / (10. * (eps + sl0mx))**2
             al1mx  = 6. / (10. * (eps + sl1mx))**2
             al2mx  = 3. / (10. * (eps + sl2mx))**2
             al0px  = 1. / (10. * (eps + sl0px))**2
             al1px  = 6. / (10. * (eps + sl1px))**2
             al2px  = 3. / (10. * (eps + sl2px))**2
             ar0mx  = 3. / (10. * (eps + sr0mx))**2
             ar1mx  = 6. / (10. * (eps + sr1mx))**2
             ar2mx  = 1. / (10. * (eps + sr2mx))**2
             ar0px  = 3. / (10. * (eps + sr0px))**2
             ar1px  = 6. / (10. * (eps + sr1px))**2
             ar2px  = 1. / (10. * (eps + sr2px))**2
             
             al0my  = 1. / (10. * (eps + sl0my))**2
             al1my  = 6. / (10. * (eps + sl1my))**2
             al2my  = 3. / (10. * (eps + sl2my))**2
             al0py  = 1. / (10. * (eps + sl0py))**2
             al1py  = 6. / (10. * (eps + sl1py))**2
             al2py  = 3. / (10. * (eps + sl2py))**2
             ar0my  = 3. / (10. * (eps + sr0my))**2
             ar1my  = 6. / (10. * (eps + sr1my))**2
             ar2my  = 1. / (10. * (eps + sr2my))**2
             ar0py  = 3. / (10. * (eps + sr0py))**2
             ar1py  = 6. / (10. * (eps + sr1py))**2
             ar2py  = 1. / (10. * (eps + sr2py))**2
             !weightings x             
             wl0mx  = al0mx / (al0mx+al1mx+al2mx)
             wl1mx  = al1mx / (al0mx+al1mx+al2mx)
             wl2mx  = al2mx / (al0mx+al1mx+al2mx)
             wl0px  = al0px / (al0px+al1px+al2px)
             wl1px  = al1px / (al0px+al1px+al2px)
             wl2px  = al2px / (al0px+al1px+al2px)
             wr0mx  = ar0mx / (ar0mx+ar1mx+ar2mx)
             wr1mx  = ar1mx / (ar0mx+ar1mx+ar2mx)
             wr2mx  = ar2mx / (ar0mx+ar1mx+ar2mx)
             wr0px  = ar0px / (ar0px+ar1px+ar2px)
             wr1px  = ar1px / (ar0px+ar1px+ar2px)
             wr2px  = ar2px / (ar0px+ar1px+ar2px)
             !weightings y
             wl0my  = al0my / (al0my+al1my+al2my)
             wl1my  = al1my / (al0my+al1my+al2my)
             wl2my  = al2my / (al0my+al1my+al2my)
             wl0py  = al0py / (al0py+al1py+al2py)
             wl1py  = al1py / (al0py+al1py+al2py)
             wl2py  = al2py / (al0py+al1py+al2py)
             wr0my  = ar0my / (ar0my+ar1my+ar2my)
             wr1my  = ar1my / (ar0my+ar1my+ar2my)
             wr2my  = ar2my / (ar0my+ar1my+ar2my)
             wr0py  = ar0py / (ar0py+ar1py+ar2py)
             wr1py  = ar1py / (ar0py+ar1py+ar2py)
             wr2py  = ar2py / (ar0py+ar1py+ar2py) 
             !negative & positive fluxes x
             flmx=vxm * (wl0mx*(c1*f(k,l,i-2,j)+c3*f(k,l,i-1,j)+c2*f(k,l,i,j))+wl1mx*(c2*f(k,l,i-1,j)+c3*f(k,l,i,j)+c1*f(k,l,i+1,j))+wl2mx*(c5*f(k,l,i,j)+c4*f(k,l,i+1,j)+c2*f(k,l,i+2,j))) 
             flpx=vxp * (wl0px*(c2*f(k,l,i-3,j)+c4*f(k,l,i-2,j)+c5*f(k,l,i-1,j))+wl1px*(c1*f(k,l,i-2,j)+c3*f(k,l,i-1,j)+c2*f(k,l,i,j))+wl2px*(c2*f(k,l,i-1,j)+c3*f(k,l,i,j)+c1*f(k,l,i+1,j)))
             frmx=vxm * (wr0mx*(c1*f(k,l,i-1,j)+c3*f(k,l,i,j)+c2*f(k,l,i+1,j))+wr1mx*(c2*f(k,l,i,j)+c3*f(k,l,i+1,j)+c1*f(k,l,i+2,j))+wr2mx*(c5*f(k,l,i+1,j)+c4*f(k,l,i+2,j)+c2*f(k,l,i+3,j)))
             frpx=vxp * (wr0px*(c2*f(k,l,i-2,j)+c4*f(k,l,i-1,j)+c5*f(k,l,i,j))+wr1px*(c1*f(k,l,i-1,j)+c3*f(k,l,i,j)+c2*f(k,l,i+1,j))+wr2px*(c2*f(k,l,i,j)+c3*f(k,l,i+1,j)+c1*f(k,l,i+2,j)))
             !negative & positive fluxes y
             flmy=vym * (wl0my*(c1*f(k,l,i,j-2)+c3*f(k,l,i,j-1)+c2*f(k,l,i,j))+wl1my*(c2*f(k,l,i,j-1)+c3*f(k,l,i,j)+c1*f(k,l,i,j+1))+wl2my*(c5*f(k,l,i,j)+c4*f(k,l,i,j+1)+c2*f(k,l,i,j+2))) 
             flpy=vyp * (wl0py*(c2*f(k,l,i,j-3)+c4*f(k,l,i,j-2)+c5*f(k,l,i,j-1))+wl1py*(c1*f(k,l,i,j-2)+c3*f(k,l,i,j-1)+c2*f(k,l,i,j))+wl2py*(c2*f(k,l,i,j-1)+c3*f(k,l,i,j)+c1*f(k,l,i,j+1)))
             frmy=vym * (wr0my*(c1*f(k,l,i,j-1)+c3*f(k,l,i,j)+c2*f(k,l,i,j+1))+wr1my*(c2*f(k,l,i,j)+c3*f(k,l,i,j+1)+c1*f(k,l,i,j+2))+wr2my*(c5*f(k,l,i,j+1)+c4*f(k,l,i,j+2)+c2*f(k,l,i,j+3)))
             frpy=vyp * (wr0py*(c2*f(k,l,i,j-2)+c4*f(k,l,i,j-1)+c5*f(k,l,i,j))+wr1py*(c1*f(k,l,i,j-1)+c3*f(k,l,i,j)+c2*f(k,l,i,j+1))+wr2py*(c2*f(k,l,i,j)+c3*f(k,l,i,j+1)+c1*f(k,l,i,j+2)))
             !fluxes x
             flx = flmx + flpx
             frx = frmx + frpx
             !fluxes y
             fly = flmy + flpy
             fry = frmy + frpy
             
             fn(i,j)  =  f(k,l,i,j) - dt/dx * (frx - flx) - dt/dy * (fry - fly)
             f(k,l,i,j) = fn(i,j)
        
        f(k,l,i,1) = fn(i,6)
        f(k,l,i,2) = fn(i,5)
        f(k,l,i,3) = fn(i,4)

        f(k,l,i,ny)  = fn(i,ny-6)
        f(k,l,i,ny-1)= fn(i,ny-5) 
        f(k,l,i,ny-2)= fn(i,ny-4)
        
        f(k,l,1,j) = fn(6,j)
        f(k,l,2,j) = fn(5,j)
        f(k,l,3,j) = fn(4,j)
        
        f(k,l,nx,j)  = fn(nx-6,j)
        f(k,l,nx-1,j)= fn(nx-5,j)
        f(k,l,nx-2,j)= fn(nx-4,j) 
    END DO
    END DO  
END DO
END DO

