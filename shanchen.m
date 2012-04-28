%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shanChen.m: Multi-component fluid, using a LB method,
%   based on the Shan-Chen model
% [X.Shan and H.Chen, http://dx.doi.org/10.1103/PhysRevE.47.1815].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, written in Matlab
% Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
% Address: EPFL-STI-LIN Station 9
% E-mail: orestis.malaspinas@epfl.ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public
% License along with this program; if not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% Boston, MA  02110-1301, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clf

% GENERAL FLOW CONSTANTS

ly           = 201;
lx           = 201;

G = -1.2;  % Amplitude of the molecular interaction force

omega1 = 1.;  % Relaxation parameter for fluid 1
omega2 = 1.;  % Relaxation parameter for fluid 2

maxT   = 80000;    % total number of iterations
tPlot  = 40;       % iterations between successive graphical outputs

% D2Q9 LATTICE CONSTANTS
tNS   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cxNS  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cyNS  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
oppNS = [  1,   4,  5,  2,  3,    8,   9,   6,   7];

[y,x] = meshgrid(1:ly,1:lx);

drho = 0.001;
delta_rho = -drho*(1-2.0*rand(lx));

% INITIAL CONDITION FOR BOTH DISTRIBUTION FUNCTIONS: (T=0) ==> TIn(i) = t(i)
for i=1:9
    fIn(i,1:lx,1:ly) = tNS(i).*(1.0 + delta_rho);
    gIn(i,1:lx,1:ly) = tNS(i).*(1.0 - delta_rho);
end

rho1  = reshape(sum(fIn),lx,ly);
imagesc(rho1');
colorbar
title('Fluid 1 density');
axis equal off; drawnow

% MAIN LOOP (TIME CYCLES)
Gomega1 = G/omega1;
Gomega2 = G/omega2;
for cycle = 1:maxT
    % MACROSCOPIC VARIABLES
    rho1 = sum(fIn);
    rho2 = sum(gIn);
    jx1  = reshape ( (cxNS * reshape(fIn,9,lx*ly)), 1,lx,ly);
    jy1  = reshape ( (cyNS * reshape(fIn,9,lx*ly)), 1,lx,ly);
    jx2  = reshape ( (cxNS * reshape(gIn,9,lx*ly)), 1,lx,ly);
    jy2  = reshape ( (cyNS * reshape(gIn,9,lx*ly)), 1,lx,ly);
   
    rhoTot_OMEGA = rho1*omega1 + rho2*omega2;
    uTotX = (jx1*omega1+jx2*omega2) ./ rhoTot_OMEGA;
    uTotY = (jy1*omega1+jy2*omega2) ./ rhoTot_OMEGA;
	
    rhoContrib1x = 0.0;
    rhoContrib2x = 0.0;
    
    rhoContrib1y = 0.0;
    rhoContrib2y = 0.0;
    for i=2:9
        rhoContrib1x = rhoContrib1x + circshift(rho1*tNS(i), [0,cxNS(i),cyNS(i)])*cxNS(i);
        rhoContrib1y = rhoContrib1y + circshift(rho1*tNS(i), [0,cxNS(i),cyNS(i)])*cyNS(i);
        
        rhoContrib2x = rhoContrib2x + circshift(rho2*tNS(i), [0,cxNS(i),cyNS(i)])*cxNS(i);
        rhoContrib2y = rhoContrib2y + circshift(rho2*tNS(i), [0,cxNS(i),cyNS(i)])*cyNS(i);
    end
    
    uTotX1 = uTotX - Gomega1.*rhoContrib2x; %POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
    uTotY1 = uTotY - Gomega1.*rhoContrib2y;
    
    uTotX2 = uTotX - Gomega2.*rhoContrib1x; %POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
    uTotY2 = uTotY - Gomega2.*rhoContrib1y;

   % COLLISION STEP FLUID 1 AND 2
   for i=1:9
      cuNS1        = 3*(cxNS(i)*uTotX1+cyNS(i)*uTotY1);
      cuNS2        = 3*(cxNS(i)*uTotX2+cyNS(i)*uTotY2);
      
      fEq(i,:,:)   = rho1 .* tNS(i) .* ...
                       ( 1 + cuNS1 + 0.5*(cuNS1.*cuNS1) - 1.5*(uTotX1.^2+uTotY1.^2) );
                       
      gEq(i,:,:)   = rho2 .* tNS(i) .* ...
                       ( 1 + cuNS2 + 0.5*(cuNS2.*cuNS2) - 1.5*(uTotX2.^2+uTotY2.^2) );
                       
      fOut(i,:,:)  = fIn(i,:,:) - omega1 .* (fIn(i,:,:)-fEq(i,:,:));
      gOut(i,:,:)  = gIn(i,:,:) - omega2 .* (gIn(i,:,:)-gEq(i,:,:));
   end

   % STREAMING STEP FLUID 1 AND 2
   for i=1:9
      fIn(i,:,:) = circshift(fOut(i,:,:), [0,cxNS(i),cyNS(i)]);
      gIn(i,:,:) = circshift(gOut(i,:,:), [0,cxNS(i),cyNS(i)]);
   end

   % VISUALIZATION
   if(mod(cycle,tPlot)==0)
       rho1     = reshape(rho1,lx,ly);
       imagesc(rho1'); colorbar
       title('Fluid 1 density');
       axis equal off; drawnow
   end
end
