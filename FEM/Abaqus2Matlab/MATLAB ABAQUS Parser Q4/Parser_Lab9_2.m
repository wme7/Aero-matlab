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
EBC1Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-1, instance=Part-1-1');
EBC2Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-2, instance=Part-1-1');
LoadS1Idx = fnFindLineWithText(lines, ...
    '*Elset, elset=_Surf-1_S1, internal, instance=Part-1-1');
LoadS2Idx = fnFindLineWithText(lines, ...
    '*Elset, elset=_Surf-1_S2, internal, instance=Part-1-1');
LoadS3Idx = fnFindLineWithText(lines, ...
    '*Elset, elset=_Surf-1_S3, internal, instance=Part-1-1');
LoadS4Idx = fnFindLineWithText(lines, ...
    '*Elset, elset=_Surf-1_S4, internal, instance=Part-1-1');
MaterialIdx = fnFindLineWithText(lines,'*Elastic');
SectionIdx = fnFindLineWithText(lines,'*Solid Section');
PressureIdx = fnFindLineWithText(lines,'*Dsload');
OutputIdx = fnFindLineWithText(lines,'*Nset, nset=Set-1, instance=Part-1-1');
SelectIdx = fnFindLineWithText(lines,'*Elset, elset=Set-1, instance=Part-1-1');

%% Import 
numNodePerElement=4;
Node=cellfun(@cellstr,lines(NodeIdx+1:ElemIdx-1));
Element=cellfun(@cellstr,lines(ElemIdx+1:EndIdx-1));
%selectOpt=str2num(lines{SelectIdx+1});

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
% Surface orientation 1
naturalBCs1=[]; surfaceOrientation1=[1];
naturalBCs1 = [naturalBCs1,str2num(lines{LoadS1Idx+1})];
naturalBCs1 = [naturalBCs1,str2num(lines{LoadS1Idx+2})];
% Surface orientation 2
naturalBCs2=[]; surfaceOrientation2=[2];
naturalBCs2 = [naturalBCs2,str2num(lines{LoadS2Idx+1})];
% Surface orientation 3
naturalBCs3=[]; surfaceOrientation3=[3];
naturalBCs3 = [naturalBCs3,str2num(lines{LoadS3Idx+1})];
% Surface orientation 4
naturalBCs4=[]; surfaceOrientation4=[4];
naturalBCs4 = [naturalBCs4,str2num(lines{LoadS4Idx+1})];

% Escential Bcs
EBC1=[]; % Nodes with fixed xDof 
EBC1=[EBC1,str2num(lines{EBC1Idx+1})];
EBC1=[EBC1,str2num(lines{EBC1Idx+2})];
EBC2=[]; % Nodes with fixed yDof 
EBC2=[EBC2,str2num(lines{EBC2Idx+1})];
EBC2=[EBC2,str2num(lines{EBC2Idx+2})];

GDof=2*numberNodes;
prescribedDof=sort([2.*EBC1-1 2.*EBC2]);
P=[0.01 0];%[Px Py]

%% Import material and section properties
Mat=str2num(lines{MaterialIdx+1});
thickness=str2num(lines{SectionIdx+1});
%DLoad=cell2str(lines{PressureIdx+1});

E=Mat(1);poisson=Mat(2);
%% Evalute force vector
forceS1=formForceVectorQ4(GDof,naturalBCs1,surfaceOrientation1,...
    elementNodes,nodeCoordinates,thickness,P);
forceS2=formForceVectorQ4(GDof,naturalBCs2,surfaceOrientation2,...
    elementNodes,nodeCoordinates,thickness,P);
forceS3=formForceVectorQ4(GDof,naturalBCs3,surfaceOrientation3,...
    elementNodes,nodeCoordinates,thickness,P);
forceS4=formForceVectorQ4(GDof,naturalBCs4,surfaceOrientation4,...
    elementNodes,nodeCoordinates,thickness,P);
force = forceS1 + forceS2 + forceS3 + forceS4;

%% Construct Stiffness matrix for T3 element
C=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);
scaleFactor=1.E6;
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

%% Compute output Edge's elements average stress,
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
% for e=1:length(outputEdge)
%     [idx idy]=find(elementNodes==outputEdge(e));
%     sigmaVM=0;
%     for k=1:length(idx)
%         temp_VMyy(k) = nodalstressyy(idx(k),idy(k));
%         temp_VMxy(k) = nodalstressxy(idx(k),idy(k));
%         temp_VMxx(k) = nodalstressxx(idx(k),idy(k));
%         sigmaVM(k)= sqrt(0.5*((temp_VMxx(k)-temp_VMyy(k))^2+...
%             (temp_VMyy(k))^2+(temp_VMxx(k))^2+6*(temp_VMxy(k))^2));
%     end
%     average_VM(e) = sum(sigmaVM)/length(sigmaVM);
% end

%% Load Abaqus results
Sxx_abaqus = importdata('Lab9_2_S11.rpt');
Uy_abaqus = importdata('Lab9_2_U2.rpt');

%% Plot Result,
y = nodeCoordinates(outputEdge,2);
disp_set = displacements(2*outputEdge);

figure(2)
subplot(1,2,1)
plot(y,disp_set,'*k');
hold on
plot(Uy_abaqus.data(:,1)+25,Uy_abaqus.data(:,2),'-b')
hold off
title('Displacement distribution');
xlabel('y (mm)');
ylabel('U_{y} (mm)');

subplot(1,2,2)
plot(y,Sxx,'*k')
hold on
plot(Sxx_abaqus.data(:,1)+25,Sxx_abaqus.data(:,2),'-b')
hold off
title('Stress distribution');
xlabel('y (mm)');
ylabel('\sigma_{xx} (MPa)');