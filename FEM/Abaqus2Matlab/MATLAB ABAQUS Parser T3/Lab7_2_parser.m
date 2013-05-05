%Finite Element Method 101-2
%National Taiwan University
%2D Elasticity problem

%%clear memory
close all; clc; clear all;
format long;

%% This function to search key word of ABAQUS inp file and returns the line number
fnFindLineWithText = @(arr,txt) ...
    find(cellfun(@(x) ~isempty (regexp(x, txt, 'once')), arr), 1);

%% Import data from Abaqus Solution,
Ux_abaqus = importdata('Ux_displacements_lab7_2.rpt');
Uy_abaqus = importdata('Uy_displacements_lab7_2.rpt');
S11_abaqus = importdata('S11_stress_lab7_2.rpt');
S22_abaqus = importdata('S22_stress_lab7_2.rpt');
S12_abaqus = importdata('S12_stress_lab7_2.rpt');
VM_abaqus = importdata('VM_stress_lab7_2.rpt');

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
EndIdx  = fnFindLineWithText(lines, '*Nset, nset=Set-1, generate');
LoadIdx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-1, instance=Part-1-1');
EBC1Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-2, instance=Part-1-1');
EBC2Idx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-3, instance=Part-1-1');
EdgeIdx = fnFindLineWithText(lines, ...
    '*Nset, nset=Set-5, instance=Part-1-1');
MaterialIdx = fnFindLineWithText(lines, '*Elastic');
SectionIdx = fnFindLineWithText(lines, '*Solid Section');
PressureIdx = fnFindLineWithText(lines, '*Dsload');

%% Import 
NodePerElement=3;
Node=cellfun(@cellstr,lines(NodeIdx+1:ElemIdx-1));
Element=cellfun(@cellstr,lines(ElemIdx+1:EndIdx-1));
numberNodes=length(Node);
numberElements=length(Element);
nodeCoordinates=zeros(numberNodes,2);
elementNodes=zeros(numberElements,NodePerElement);

%% Import nodes and elements
for i=1:numberNodes
    temp=str2num(Node{i});
    nodeCoordinates(i,:)=temp(2:end);
end

for i=1:numberElements
    temp=str2num(Element{i});
    elementNodes(i,:)=temp(2:end);
end

drawingMesh(nodeCoordinates,elementNodes,'T3','b-o');

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

%% Import Edge nodes
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
%force=formForceVectorT3(GDof,naturalBCs,surfaceOrientation,...
%    elementNodes,nodeCoordinates,P,thickness);
force = zeros(GDof,1);
force(loadDof) = P;

%% Construct Stiffness matrix for T3 element
C=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% matrix D
D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

%% stress
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    elementDof = [];
    for i = elementNodes(e,:)
        temp = [2*i-1,2*i];
        elementDof = [elementDof,temp];
    end
    B = Bmat(nodeCoordinates,elementNodes,e);
    stress(:,e) = D*B*displacements(elementDof);
end

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);

scaleFactor=1.E3;
drawingMesh(nodeCoordinates+scaleFactor*[displacements(1:2:2*numberNodes) ...
    displacements(2:2:2*numberNodes)],elementNodes,'T3','r--');

%% Find elements arround every nodes in the vertical edge: EdgeNodes,
for i = 1:length(EdgeNodes)
    [row,col] = find(elementNodes==EdgeNodes(i));
    elementList{i} = row;
end

%% Compute stress averages for the nodes in EBC1
stress = stress';
for i = 1:length(EdgeNodes)
    s11 = stress(elementList{i},1);
    s11_bar(i) = mean(s11);
    s22 = stress(elementList{i},2);
    s22_bar(i) = mean(s22);
    s12 = stress(elementList{i},3);
    s12_bar(i) = mean(s12);
    Se = sqrt(s11.^2 + s22.^2 -s11.*s22 + 3*s12.*s12); %Von Mises stress
    Se_bar(i) = mean(Se);
end

%% Compute local coodinate y-Position of target Edge nodes,
y = nodeCoordinates(EdgeNodes,2);

%% Compute y-displacements in EBC1 nodes
Ux = displacements(2*EdgeNodes-1);
Uy = displacements(2*EdgeNodes);

%% plot data
figure(2)

% Displacements Ux,
subplot(3,2,1); hold on;
plot(y,Ux,'*b'); 
plot(sort(y),Ux_abaqus.data(:,2),'-r'); 
title('Displacements'); xlabel('U_1,(mm)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',4);
hold off;

% Displacements Uy,
subplot(3,2,2); hold on;
plot(y,Uy,'*b'); 
plot(sort(y),Uy_abaqus.data(:,2),'-r'); 
title('Displacements'); xlabel('U_2,(mm)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',1);
hold off;

% Stress S11,
subplot(3,2,3); hold on; 
plot(y,s11_bar,'*b'); 
plot(sort(y),S11_abaqus.data(:,2),'-r');
title('Stress distribution'); xlabel('\sigma_11,(MPa)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',2);
hold off

% Stress S22,
subplot(3,2,4); hold on; 
plot(y,s22_bar,'*b'); 
plot(sort(y),S22_abaqus.data(:,2),'-r');
title('Stress distribution'); xlabel('\sigma_22,(MPa)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',2);
hold off

% Stress S12,
subplot(3,2,5); hold on; 
plot(y,s12_bar,'*b'); 
plot(sort(y),S12_abaqus.data(:,2),'-r'); 
title('Stress distribution'); xlabel('\sigma_12,(MPa)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',4);
hold off

% Stress VM,
subplot(3,2,6); hold on; 
plot(y,Se_bar,'*b'); 
plot(sort(y),VM_abaqus.data(:,2),'-r'); 
title('Stress distribution'); xlabel('\sigma_Mises,(MPa)'); ylabel('y,(mm)');
legend('Matlab','Abaqus',1);
hold off