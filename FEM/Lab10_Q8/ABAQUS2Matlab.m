%Finite Element Method 101-2
%National Taiwan University
%2D Elasticity problem

%%clear memory
close all; clc; clear all;
%% This function to search key word of ABAQUS inp file and returns the line number
% fnFindLineWithText = @(arr,txt) ...
%     find(cellfun(@(x) ~isempty(regexp(x,txt, 'once')), arr), 1);
% %% Open file with user interface
% [filename,filepath]=uigetfile('*.inp','Select Input file');
% file = [filepath filename];
% fileID = fopen(file,'r');
% if fileID == -1
%     disp('Error opening the file')
% end
% 
% %% Import all lines into a cell
% lines = {};
% while ~feof(fileID)
%     line = fgetl(fileID);
%     lines{end+1} = line;
% end
% 
% fclose(fileID);
% 
% %% this part to locate key word
% NodeIdx = fnFindLineWithText(lines, '*Node');
% ElemIdx = fnFindLineWithText(lines, '*Element');
% EndIdx = fnFindLineWithText(lines, '*Nset, nset=Set-1, generate');
% EBCIdx = fnFindLineWithText(lines, ...
%     '*Nset, nset=Set-1, instance=Part-1-1, generate');
% MaterialIdx = fnFindLineWithText(lines,'*Elastic');
% SectionIdx = fnFindLineWithText(lines,'*Solid Section');
% PressureIdx = fnFindLineWithText(lines,'*Dsload');
% OutputIdx = fnFindLineWithText(lines,'*Nset, nset=output, generate');
% SelectIdx = fnFindLineWithText(lines,'*Elset, elset=output');
% 
% %% Import 
% numNodePerElement=4;
% Node=cellfun(@cellstr,lines(NodeIdx+1:ElemIdx-1));
% Element=cellfun(@cellstr,lines(ElemIdx+1:EndIdx-1));
% selectOpt=str2num(lines{SelectIdx+1});
% 
% numberNodes=length(Node);
% numberElements=length(Element);
% nodeCoordinates=zeros(numberNodes,2);
% elementNodes=zeros(numberElements,numNodePerElement);
% %% Import nodes and elements
% for i=1:numberNodes
%     temp=str2num(Node{i});
%     nodeCoordinates(i,:)=temp(2:end);
% end
% 
% for i=1:numberElements
%     temp=str2num(Element{i});
%     elementNodes(i,:)=temp(2:end);
% end
% 
% drawingMesh(nodeCoordinates,elementNodes,'Q4','b-o');
%% Load Mesh for Lab 10
[nodeCoordinates,elementNodes]=Mesh_Lab10('Q4');
% node coordinates are given in mm
NodePerElement=4;
numberNodes=size(nodeCoordinates,1);
numberElements=size(elementNodes,1);

%% Import BCs
% naturalBCs=[];
% surfaceOrientation=[];
% 
% for i=1:numNodePerElement
%     LoadIdx = fnFindLineWithText(lines,...
%         ['**Elset, elset=_Surf-1_S',num2str(i)]);
%     if ~isempty(LoadIdx)
%         naturalBCs=[naturalBCs;str2num(lines{LoadIdx+1})];
%         surfaceOrientation=[surfaceOrientation;i];
%     end
% end
% 
% generateBC=str2num(lines{EBCIdx+1});
% EBC=generateBC(1):generateBC(3):generateBC(2);
% outTemp_1=str2num(lines{OutputIdx+1});
% outputEdge=outTemp_1(1):outTemp_1(3):outTemp_1(2);
% GDof=2*numberNodes;
% prescribedDof=sort([2.*EBC-1 2.*EBC]);
% P=[0 -20];%[Px Py]

% Essential BC's
GDof = 2*numberNodes;
prescribedDof = [1,2,3,4];

% Natural BC's
force = zeros(GDof,1);
force(end) = -10; % 10 [kN]

%% Import material and section properties
% Mat=str2num(lines{MaterialIdx+1});
% thickness=str2num(lines{SectionIdx+1});
% %DLoad=cell2str(lines{PressureIdx+1});
% 
% E=Mat(1);poisson=Mat(2);

E = 3E7; % [GPa]
poisson = 0.3; %[-]
thickness = 1; % [mm]
%% Evalute force vector
%force=formForceVectorQ4(GDof,naturalBCs,surfaceOrientation,...
%    elementNodes,nodeCoordinates,thickness);
%% Construct Stiffness matrix for Q8 element
C=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];

stiffness=formStiffness2D(GDof,numberElements,...
    elementNodes,numberNodes,nodeCoordinates,C,thickness);

%% solution
displacements=solution(GDof,prescribedDof,stiffness,force);

%% output displacements
outputDisplacements(displacements, numberNodes, GDof);
scaleFactor=1.E8;
drawingMesh(nodeCoordinates+scaleFactor*[displacements(1:2:2*numberNodes) ...
    displacements(2:2:2*numberNodes)],elementNodes,'Q4','r--');
