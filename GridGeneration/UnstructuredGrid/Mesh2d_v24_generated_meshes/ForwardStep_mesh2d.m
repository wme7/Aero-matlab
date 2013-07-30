%*************************************************************************
% Forward step mesh subroutine using mesh2d with
% -> use drawingMesh.p to verify
% coded by Manuel Diaz, NTU, 2013.07.30
%*************************************************************************

% Forward Step Domain Triangular Mesh

clear all; clc;

% file name
name = 'ForwardStep.plt';

%% Define Reference nodes and their connectivity Matrix
nodes = [0.0 0.0;
        0.6 0.0;
        0.6 0.2;
        3.0 0.2;
        3.0 1.0;
        0.0 1.0];
    
cnect = [1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 1];

a = 0.6; b = 0.2; 
h = @(x,y) 1/35 + 1/25*sqrt( (x-a).^2+(y-b).^2 );
    
% Element Size, hdata
hdata = [];                 % Initialize hdata class
hdata.hmax = 1/25;          % Domain element size.
%hdata.edgeh = [2,1/20];     % Element size on specified geometry edges.
hdata.fun = h;              % User defined size function.

% Solver options
%options.dhmax = 0.01;

% Build mesh
[p,t] = mesh2d(nodes,cnect,hdata);
[p,t] = fixmesh(p,t);

%% Mesh data for DGFEM
K=size(t,1); % Number of elements
Nv=size(p,1); % Number of vertices in mesh
Nfaces=size(t,2); % Number of faces/element
VX = p(:,1); % Vertice x-coordinates
VY = p(:,2); % Vertice y-coordinates
EToV = t; % Element to Vertice table

% plot with drawingmesh.p
drawingMesh([VX VY],EToV,'T3','k-o');

% Clear extra data
clear('nodes','cnect','t','a','b','h','p','hdata');

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
