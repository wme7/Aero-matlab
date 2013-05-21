%Finite Element Method 101-2
%National Taiwan University
%2D Elasticity problem

%%clear memory
close all; clc; clear all;
%% This function to search key word of ABAQUS inp file and returns the line number
fnFindLineWithText = @(arr,txt) ...
    find(cellfun(@(x) ~isempty(regexp(x,txt, 'once')), arr), 1);
%% Open file with user interface
[filename,filepath]=uigetfile('*.inp','Select Input file');
file = [filepath filename];
fileID = fopen(file,'r');
if fileID == -1
    disp('Error opening the file')
end

%% Import all lines into a cell
lines = {};
while ~feof(fileID)
    line = fgetl(fileID);
    lines{end+1} = line;
end

fclose(fileID);

%% this part to locate key word
NodeIdx = fnFindLineWithText(lines, '*Node');
ElemIdx = fnFindLineWithText(lines, '*Element');
EndIdx = fnFindLineWithText(lines, '*Nset, nset=Set-1, generate');
LoadIdx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-1, instance=Part-1-1');
EBC1Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-2, instance=Part-1-1');
EBC2Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-3, instance=Part-1-1');
EdgeIdx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-4, instance=Part-1-1');
MaterialIdx = fnFindLineWithText(lines,'*Elastic');
SectionIdx = fnFindLineWithText(lines,'*Solid Section');
PressureIdx = fnFindLineWithText(lines,'*Dsload');
OutputIdx = fnFindLineWithText(lines,'*Nset, nset=Set-4, instance=Part-1-1');

%% Import 
numNodePerElement=4;
Node=cellfun(@cellstr,lines(NodeIdx+1:ElemIdx-1));
Element=cellfun(@cellstr,lines(ElemIdx+1:EndIdx-1));

numberNodes=length(Node);
numberElements=length(Element);
nodeCoordinates=zeros(numberNodes,2);
elementNodes=zeros(numberElements,numNodePerElement);

%% Import nodes and elements
for idy=1:numberNodes
    temp=str2num(Node{idy});
    nodeCoordinates(idy,:)=temp(2:end);
end

for idy=1:numberElements
    temp=str2num(Element{idy});
    elementNodes(idy,:)=temp(2:end);
end

drawingMesh(nodeCoordinates,elementNodes,'Q4','b-o');

%% Import BCs
naturalBCs=[]; % node of the concentrated force
surfaceOrientation=[];

generateBC=cellfun(@cellstr,lines(LoadIdx+1));
for i = 1:length(generateBC)
    temp=str2num(generateBC{i});
    naturalBCs=[naturalBCs,temp];
end

EBC1 = []; % fixed displacements nodes in x and y

generateBC=cellfun(@cellstr,lines(EBC1Idx+1));
for i = 1:length(generateBC)
    temp=str2num(generateBC{i});
    EBC1=[EBC1,temp];
end

EBC2 = []; % fixed displacements nodes in x and y

generateBC=cellfun(@cellstr,lines(EBC2Idx+1));
for i = 1:length(generateBC)
    temp=str2num(generateBC{i});
    EBC2=[EBC2,temp];
end

EBC = [EBC1,EBC2];

GDof=2*numberNodes;
prescribedDof=sort([2.*EBC-1 2.*EBC]);
loadDof=sort([2.*naturalBCs-1 2.*naturalBCs]);
P=[700 0];%[Px Py]

%% Import output Edge nodes
EdgeNodes = []; 

generateBC=cellfun(@cellstr,lines(EdgeIdx+1));
for i = 1:length(generateBC)
    temp=str2num(generateBC{i});
    EdgeNodes=[EdgeNodes,temp];
end

%% Import material and section properties
Mat=str2num(lines{MaterialIdx+1});
thickness=str2num(lines{SectionIdx+1});
%DLoad=cell2str(lines{PressureIdx+1});

E=Mat(1);poisson=Mat(2);

%% Evalute force vector
force = zeros(GDof,1);
force(loadDof) = P;

%% Construct Stiffness matrix for T3 element
C=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);
scaleFactor=1.E3;
drawingMesh(nodeCoordinates+scaleFactor*[displacements(1:2:2*numberNodes) ...
    displacements(2:2:2*numberNodes)],elementNodes,'Q4','r--');

%% Using plane-stress 
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

%% Displacement of Output nodes
outputEdge =[]; 
outputEdge = [outputEdge,str2num(lines{OutputIdx+1})];
outputEdge = [outputEdge,str2num(lines{OutputIdx+2})];

%% Compute B Matrix and Strain for each element.
% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gauss2d('2x2');

for e=1:numberElements
    numNodePerElement = length(elementNodes(e,:));
    numEDOF = 2*numNodePerElement;
    elementDof=zeros(1,numEDOF);
    for i = 1:numNodePerElement
        elementDof(2*i-1)=2*elementNodes(e,i)-1;
        elementDof(2*i)=2*elementNodes(e,i);
    end
  
    % cycle for Gauss point
    for q=1:size(gaussWeights,1)     % for 2x2 case, q: 1,2,3,4
        GaussPoint=gaussLocations(q,:);
        xi=GaussPoint(1);
        eta=GaussPoint(2);

        % shape functions and derivatives
        [shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);

        % Jacobian matrix, inverse of Jacobian,
        % derivatives w.r.t. x,y
        [Jacob,invJacobian,XYderivatives]=...
            Jacobian(nodeCoordinates(elementNodes(e,:),:),naturalDerivatives);

        %  B matrix for element 'e'
        B=zeros(3,numEDOF);
        B(1,1:2:numEDOF)  =     XYderivatives(:,1)';
        B(2,2:2:numEDOF)  = XYderivatives(:,2)';
        B(3,1:2:numEDOF)  =     XYderivatives(:,2)';
        B(3,2:2:numEDOF)  = XYderivatives(:,1)';

        %  element 'e' displacement vector (for Q4: 8-Dofs)
        d(1:2:8) = displacements([2*elementNodes(e,:)-1],1); % x-disp
        d(2:2:8) = displacements([2*elementNodes(e,:)],1);   % y-disp

        %  element strain & stress vectors
        strain = B*d';
        stress(:,q) = D*strain;
    end  
    
    % Another cycle for Gauss points
    % Computing shape functions using (r,s) coords,
    for q=1:size(gaussWeights,1)     % for 2x2 case, q: 1,2,3,4
        GaussPoint=gaussLocations(q,:);
        r=1/GaussPoint(1);
        s=1/GaussPoint(2);
        
    % shape functions and derivatives
        [shapeFunction,naturalDerivatives]=shapeFunctionQ4(r,s);
    
    % interpolate stress using new shape functions to nodal points,
    % S_P(r,s) = N(r,s)*S_GP
    % Sxx(r,s) = N1*S1xx + N2*S2xx + N3*S3xx + N4*S4xx
    % Syy(r,s) = N1*S1yy + N2*S2yy + N3*S3yy + N4*S4yy
    % Sxy(r,s) = N1*S1xy + N2*S2xy + N3*S3xy + N4*S4xy
        nodalstressxx(e,q) = stress(1,:)*shapeFunction;
        nodalstressyy(e,q) = stress(2,:)*shapeFunction;
        nodalstressxy(e,q) = stress(3,:)*shapeFunction;
    
    end     
end

%% Compute output Edge's nodes average stress,
for e = 1:length(outputEdge)
    [idx idy] = find(elementNodes==outputEdge(e));
    anSyy=0;     anSxy=0;     anSxx=0;
    for k=1:length(idx)
        temp_Syy = nodalstressyy(idx(k),idy(k));
        anSyy = temp_Syy + anSyy;
        temp_Sxy = nodalstressxy(idx(k),idy(k));
        anSxy = temp_Sxy + anSxy;
        temp_Sxx = nodalstressxx(idx(k),idy(k));
        anSxx = temp_Sxx + anSxx;
    end
    Syy(e)= anSyy/length(idx);
    Sxy(e)= anSxy/length(idx);
    Sxx(e)= anSxx/length(idx);
end

%% Von Mises
for e=1:length(outputEdge)
    [idx idy]=find(elementNodes==outputEdge(e));
    sigmaVM=0;
    for k=1:length(idx)
        temp_VMyy(k) = nodalstressyy(idx(k),idy(k));
        temp_VMxy(k) = nodalstressxy(idx(k),idy(k));
        temp_VMxx(k) = nodalstressxx(idx(k),idy(k));
        sigmaVM(k)= sqrt(0.5*((temp_VMxx(k)-temp_VMyy(k))^2+...
            (temp_VMyy(k))^2+(temp_VMxx(k))^2+6*(temp_VMxy(k))^2));
    end
    average_VM(e) = sum(sigmaVM)/length(sigmaVM);
end

%% Load Abaqus results
Sxx_abaqus = importdata('Lab9_3_S11.rpt');
Sxy_abaqus = importdata('Lab9_3_S12.rpt');
Syy_abaqus = importdata('Lab9_3_S22.rpt');
Ux_abaqus  = importdata('Lab9_3_U1.rpt');
Uy_abaqus  = importdata('Lab9_3_U2.rpt');
VM_abaqus  = importdata('Lab9_3_VM.rpt');

%% Plot Result,
y = nodeCoordinates(outputEdge,2);
disp_x = displacements(2*outputEdge-1);
disp_y = displacements(2*outputEdge);

figure(2)
subplot(2,3,1)
plot(y,disp_x,'*k');
title('Displacement distribution');
hold on
plot(Ux_abaqus.data(:,1),Ux_abaqus.data(:,2),'-b')
hold off
xlabel('y (mm)');
ylabel('U_{x} (mm)');

subplot(2,3,2)
plot(y,disp_y,'*k');
title('Displacement distribution');
hold on
plot(Uy_abaqus.data(:,1),Uy_abaqus.data(:,2),'-b')
hold off
xlabel('y (mm)');
ylabel('U_{y} (mm)');

subplot(2,3,3)
plot(y,Sxx,'*k')
hold on
plot(Sxx_abaqus.data(:,1),Sxx_abaqus.data(:,2),'-b')
hold off
title('Stress distribution');
xlabel('y (mm)');
ylabel('\sigma_{xx} (MPa)');

subplot(2,3,4)
plot(y,Syy,'*k')
hold on
plot(Syy_abaqus.data(:,1),Syy_abaqus.data(:,2),'-b')
hold off
title('Stress distribution');
xlabel('y (mm)');
ylabel('\sigma_{yy} (MPa)');

subplot(2,3,5)
plot(y,Sxy,'*k')
hold on
plot(Sxy_abaqus.data(:,1),Sxy_abaqus.data(:,2),'-b')
hold off
title('Stress distribution');
xlabel('y (mm)');
ylabel('\sigma_{xy} (MPa)');

subplot(2,3,6)
plot(y,average_VM,'*k')
hold on
plot(VM_abaqus.data(:,1),VM_abaqus.data(:,2),'-b')
hold off
title('Stress distribution');
xlabel('y (mm)');
ylabel('\sigma_Mises (MPa)');