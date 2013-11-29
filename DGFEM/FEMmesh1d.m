%% 1D grid
% I'll be following ideas from Allan P. Engsig-Karup lecture notes.
clear all; clc;

% Let's build an equidistant 1d grid of elements
a=0; b=1; nE=10; nV=nE+1; nF=2; nFp=1; pDeg=5; eNodes=pDeg+1; nodetol=1e-6;

% Vertex Coordinates VX, VY, VZ
VX = linspace(a,b,nV);

% Element to Vertices: EToV
EtoV = zeros(nE,nF);
for e = 1:nE
    EtoV(e,:) = [e,e+1];
end

% Nodes Coordiantes of Reference element
r = [-1,-1/4,-1/2,1/2,1/3,1];

% build x grid nodes
va = EtoV(:,1); vb = EtoV(:,2);
dx = VX(vb)-VX(va); 
xc = (VX(vb)+VX(va))/2;
x = ones(size(r))'*xc + r'*dx/2;

% Jacobian
J = (x(:,2:end)-x(:,1:end-1))/2;

% Compute outward pointing normals at element faces
nx = zeros(nFp*nF,nE); %nFp number of Nodes on Face
nx(1,:)=-1; nx(2,:)= 1;

% Correct Table EToV
    %no needed in 1D

% total faces
TotalFaces = nF*nE;

% local face to local node connections
FtoN = [1,2];

FtoV = spalloc(TotalFaces,nV,2*TotalFaces);
se=1; % element surface
for e = 1:nE
    for f = 1:nF
        FtoV(se,EtoV(e,FtoN(f))) = 1;
        se = se + 1;
    end
end
%spy(FtoV);

% face to face connectivity
FtoF = FtoV*FtoV' - speye(TotalFaces); spy(FtoF)
[faces1,faces2] = find(FtoF==1);

% faces to element and face numbers
element1 = floor( (faces1-1)/nF )+1;
face1    =   mod( (faces1-1),nF )+1;
element2 = floor( (faces2-1)/nF )+1;
face2    =   mod( (faces2-1),nF )+1;

% Element to faces & Element to Element
% this arrays contain information of the neighbor element and faces
EtoE = (1:nE)'*ones(1,nF);
EtoF = ones(1,nE)'*(1:nF);
ind = sub2ind(size(EtoF),element1,face1);
EtoE(ind) = element2;
EtoF(ind) = face2;

% build node ids
nodeids = reshape((1:nE*eNodes),eNodes,nE);

% Access element at the faces of the reference element
fmask1 = find(abs(r+1) < nodetol);
fmask2 = find(abs(r-1) < nodetol);
Fmask = [fmask1;fmask2]';
Fx = x(Fmask(:),:); % Face x nodes
FinvJ = 1./J(Fmask,:);% Face Jacobian

% Build Maps
vmapM = zeros(nFp,nF,nE);
vmapP = zeros(nFp,nF,nE);
for e1 = 1:nE
    for f1 = 1:nF
        % find index of face nodes with respect to volume node ordering
        vmapM(:,f1,e1) = nodeids(Fmask(:,f1),e1);
    end
end

for e1 = 1:nE
    for f1 = 1:nF
        % find neighbor
        e2 = EtoE(e1,f1); f2 = EtoF(e1,f1);
        
        % find volume node numbers of the left(M) and right(P) nodes
        vidM = vmapM(:,f1,e1); vidP = vmapM(:,f2,e2);
        x1 = x(vidM); x2 = x(vidP);
        
        % Compute distance matrix
        Dist = (x1-x2).^2;
        if Dist < nodetol; vmapP(:,f1,e1) = vidP; end
    end
end
% reduce the size of maps
vmapP = vmapP(:); vmapM = vmapM(:);

% list of boundary nodes
mapB = find(vmapP==vmapM); vmapB = vmapP(mapB);

% create specific left (inflow) and right (outflow) maps
mapI = 1; mapO = nE*nF; vmapI = vmapM(mapI); vmapO = vmapP(mapO);
