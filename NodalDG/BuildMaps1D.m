function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D

% function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D
% Purpose: Connectivity and boundary tables for nodes given in the Nv # of elements,
% 	       each with N+1 degrees of freedom.

Globals1D;

% number volume nodes consecutively
nodeids = reshape(1:Nv*Np, Np, Nv);
vmapM   = zeros(Nfp, Nfaces, Nv); 
vmapP   = zeros(Nfp, Nfaces, Nv); 

for k1=1:Nv
  for f1=1:Nfaces
    % find index of face nodes with respect to volume node ordering
    vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
  end
end

for k1=1:Nv
  for f1=1:Nfaces
    % find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
    
    % find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
    
    x1  = x(vidM); x2  = x(vidP);
    
    % Compute distance matrix
    D = (x1 -x2 ).^2;
    if (D<NODETOL) vmapP(:,f1,k1) = vidP; end;
  end
end

vmapP = vmapP(:); vmapM = vmapM(:);

% Create list of boundary nodes
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);

% Create specific left (inflow) and right (outflow) maps
mapI = 1; mapO = Nv*Nfaces; vmapI = 1; vmapO = Nv*Np;
return
