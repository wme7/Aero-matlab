%*************************************************************************
% Forward step mesh subroutine using mesh2d with
% -> use drawingMesh.p to verify
% coded by Manuel Diaz, NTU, 2013.07.30
%*************************************************************************

% Forward Step Domain Triangular Mesh

%%clear memory
close all; clc; clear all; format long;

% file name
name = 'ForwardStep.plt';

%% This function to search key word of ABAQUS inp file and returns the line number
fnFindLineWithText = @(arr,txt) ...
    find(cellfun(@(x) ~isempty (regexp(x, txt, 'once')), arr), 1);

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
EBCIdx = fnFindLineWithText(lines, ...
    '*Nset, nset=EBC, instance=Plate-1, generate')
MaterialIdx = fnFindLineWithText(lines,'*Elastic');
SectionIdx = fnFindLineWithText(lines,'*Solid Section');
PressureIdx = fnFindLineWithText(lines,'*Dsload');

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

%% Capture Data for DGFEM
K = numberElements; % Number of elements
Nv = numberNodes; % Number of vertices in mesh
Nfaces=size(elementNodes,2); % Number of faces/element
VX = nodeCoordinates(:,1); % Vertice x-coordinates
VY = nodeCoordinates(:,2); % Vertice y-coordinates
EToV = elementNodes; % Element to Vertice table

% plot with drawingmesh.p
drawingMesh(nodeCoordinates,elementNodes,'T3','b-o');

% create fictitious data
Z = ones(size(VX));

%% Export to Tecplot
file = fopen(name,'w');
% h1 gets the handel for the file "mesh1.tec".
% 'w' specifies that it will be written.
% similarly 'r' is for reading and 'a' for appending.
fprintf(file, 'TITLE = Forward Step Mesh"\n');
fprintf(file, 'VARIABLES = "X", "Y", "Z"\n');

fprintf(file, 'ZONE n=%d, e=%d, f=fepoint, et=triangle \n\n',Nv,K);
for i = 1:Nv % for all vertices
	fprintf(file, '%f\t%f\t%f\n', VX(i,1),VY(i,1),Z(i));
end
fprintf(file, '\n'); % space
for i = 1:K % for all elements
    fprintf(file, '%i\t',EToV(i,:)); fprintf(file,'\n');
end

fclose(file);