%% 2D grid
clear all; close all; clc;

% Select mesh case
meshcase = 'Distmesh';

switch meshcase
    case 'Q4' % Quadrilateral mesh
        Ex=8; Ey=4; 
        [VX VY EtoV nV nE] = Q4meshGenerator([0,1],[0,1],Ex,Ey);
    case 'T3' % Triangular mesh
        Ex=8; Ey=4; 
        [VX VY EtoV nV nE] = T3meshGenerator([0,1],[0,1],Ex,Ey);
    case 'AG' % Analytic Grids
        Ex=5; Ey=10; 
        [VX VY EtoV nV nE] = AGmeshGenerator([0,1],[0,1],Ex,Ey,'EllipticCylinderCoordinates');
        %[VX VY EtoV nV nE] = AGmeshGenerator([0,1],[0,1],Ex,Ey,'PolarCoordinates');
    case 'Distmesh' % DistMesh
        [VX VY EtoV nV nE] = DistMeshGenerator();
    case 'mesh2d' % Mesh2D
        [VX VY EtoV nV nE] = M2DmeshGenerator();
end

% Correct Table EToV
switch meshcase
    case {'T3','Distmesh','mesh2d'}
        % Reorder elements to ensure counter clockwise orientation
        ax = VX(EtoV(:,1)); ay = VY(EtoV(:,1));
        bx = VX(EtoV(:,2)); by = VY(EtoV(:,2));
        cx = VX(EtoV(:,3)); cy = VY(EtoV(:,3));
        D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
        i = find(D<0);
        EtoV(i,:) = EtoV(i,[1 3 2]);
    otherwise
        % Correct elements
        ax = VX(EtoV(:,1)); ay = VY(EtoV(:,1));
        bx = VX(EtoV(:,2)); by = VY(EtoV(:,2));
        cx = VX(EtoV(:,3)); cy = VY(EtoV(:,3));
        dx = VX(EtoV(:,4)); dy = VY(EtoV(:,4));
        Area = 0.5*(-ay.*bx+ax.*by-by.*cx+bx.*cy-cy.*dx+cx.*dy);
        %Area = 0.5*(ax-dx).*(cy-dy)-(bx-dx).*(by-dy)-(cx-dx).*(ay-dy);
        i = find(Area<0); fprintf('problem elements: %1.0f\n',i);
        EtoV(i,:) = EtoV(i,[1 4 3 2]);
end

% Visualize elements
figure(2)
    patch('Faces',EtoV,'Vertices',[VX,VY],...
        'FaceColor','none','EdgeColor','k'); hold on;   
    
        % Plot nodes and elements info
        for i = 1:nV
            text(VX(i),VY(i),int2str(i),'fontsize',8,....
                'fontweight','bold','Color','r'); 
        end
        for e = 1:nE
            pos = [mean(VX(EtoV(e,:))),mean(VY(EtoV(e,:)))];
            text(pos(1),pos(2),int2str(e),'fontsize',8,...
                'fontweight','bold','color','b');
        end

% Basic information
nF=2; nFp=1; pDeg=5; eNodes=pDeg+1; nodetol=1e-6;

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
