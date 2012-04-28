%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cavity2d.m: 2D cavity flow, simulated by a LB method            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, Matlab script
% Copyright (C) 2006-2008 Jonas Latt
% Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
% E-mail: Jonas.Latt@cui.unige.ch
%
% Implementation of 2d cavity geometry and Zou/He boundary
% condition by Adriano Sciacovelli
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

% GENERAL FLOW CONSTANTS 
lx = 128; 
ly = 128; 

uLid  = 0.05; % horizontal lid velocity 
vLid  = 0;    % vertical lid velocity 
Re    = 100;  % Reynolds number 
nu    = uLid *lx / Re;     % kinematic viscosity 
omega = 1. / (3*nu+1./2.); % relaxation parameter 
maxT  = 40000; % total number of iterations 
tPlot = 10;    % cycles for graphical output

% D2Q9 LATTICE CONSTANTS 
t   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; 
cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1]; 
opp = [ 1, 4, 5, 2,  3, 8, 9,  6,  7]; 
lid = [2: (lx-1)]; 

[y,x] = meshgrid(1:ly,1:lx); 
obst = ones(lx,ly); 
obst(lid,2:ly) = 0; 
bbRegion = find(obst); 

% INITIAL CONDITION: (rho=0, u=0) ==> fIn(i) = t(i) 
fIn = reshape( t' * ones(1,lx*ly), 9, lx, ly); 

% MAIN LOOP (TIME CYCLES) 
for cycle = 1:maxT 

% MACROSCOPIC VARIABLES 
rho = sum(fIn); 
ux = reshape ( (cx * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 
uy = reshape ( (cy * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 

% MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS 

ux(:,lid,ly) = uLid; %lid x - velocity 
uy(:,lid,ly) = vLid; %lid y - velocity 
rho(:,lid,ly) = 1 ./ (1+uy(:,lid,ly)) .* ( ... 
                sum(fIn([1,2,4],lid,ly)) + 2*sum(fIn([3,6,7],lid,ly)) ); 

% MICROSCOPIC BOUNDARY CONDITIONS: LID (Zou/He BC)
fIn(5,lid,ly) = fIn(3,lid,ly) - 2/3*rho(:,lid,ly).*uy(:,lid,ly); 
fIn(9,lid,ly) = fIn(7,lid,ly) + 1/2*(fIn(4,lid,ly)-fIn(2,lid,ly))+ ... 
                1/2*rho(:,lid,ly).*ux(:,lid,ly) - 1/6*rho(:,lid,ly).*uy(:,lid,ly); 
fIn(8,lid,ly) = fIn(6,lid,ly) + 1/2*(fIn(2,lid,ly)-fIn(4,lid,ly))- ... 
                1/2*rho(:,lid,ly).*ux(:,lid,ly) - 1/6*rho(:,lid,ly).*uy(:,lid,ly); 

% COLLISION STEP 
for i=1:9 
    cu = 3*(cx(i)*ux+cy(i)*uy); 
    fEq(i,:,:) = rho .* t(i) .* ... 
        ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) ); 
    fOut(i,:,:) = fIn(i,:,:) - omega .* (fIn(i,:,:)-fEq(i,:,:)); 
end 

% MICROSCOPIC BOUNDARY CONDITIONS: NO-SLIP WALLS (bounce-back)
for i=1:9 
    fOut(i,bbRegion) = fIn(opp(i),bbRegion); 
end 

% STREAMING STEP 
for i=1:9 
    fIn(i,:,: ) = circshift(fOut(i,:,: ), [0,cx(i),cy(i)]); 
end

% VISUALIZATION
if (mod(cycle,tPlot)==0)
    u = reshape(sqrt(ux.^2+uy.^2),lx,ly);
    u(bbRegion) = nan;
    imagesc(u(:,ly:-1:1)'./uLid);
    colorbar
    axis equal off; drawnow
end

end 
