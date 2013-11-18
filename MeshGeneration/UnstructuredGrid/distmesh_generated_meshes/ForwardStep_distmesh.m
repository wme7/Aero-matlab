%*************************************************************************
% Forward step mesh subroutine using distmesh 
% -> fixmesh subroutine is implemented to merge the repeated nodes
% coded by Manuel Diaz, NTU, 2013.07.30
%*************************************************************************

% Forward Step Domain Triangular Mesh

clear all; clc;

% file name
name = 'ForwardStep.plt';

%% First Method: using dpoly
% pv = [0 0;0.6 0;0.6 0.2; 3 0.2;3 1;0 1];
% [p,t] = distmesh2d(@dpoly,@huniform,0.125,[0,0; 3,1],pv,pv); %FAIL!

%% Second Method: using the difference of two shapes
fd=inline('ddiff(drectangle(p,0,3,0,1),drectangle(p,0.6,3,0,0.2))','p');
pfix = [0 0;0.6 0;0.6 0.2; 3 0.2;3 1;0 1];
    % Build Mesh
    [p,t] = distmesh2d(fd,@huniform,1/20,[0,0; 3,1],pfix);
    % Fix output data
    [p,t] = fixmesh(p,t);

%% Mesh data for DGFEM
K=size(t,1); % Number of elements
Nv=size(p,1); % Number of vertices in mesh
Nfaces=size(t,2); % Number of faces/element
VX = p(:,1); % Vertice x-coordinates
VY = p(:,2); % Vertice y-coordinates
EToV = t; % Element to Vertice table

% Clear extra data
clear('pfix','t');

%% Visualization
close all;
% Use drawingMesh (FEM) function
%figure(1); drawingMesh([VX VY],EToV,'T3','k-o');

% Use dist_plot to visualize the distance function:
%figure(2); dist_plot([VX VY],EToV,fd);

% Using Matlab Triangular Plot
triplot(EToV,VX,VY,'k')

% %% Create a List of boundary nodes
% fdInner = inline(dcircle(p,0,0,0.4),p);
% nodesInner = find(abs(fdInner([p]))<1e-3);
% fdOuter = inline(drectangle(p,-1,1,-1,1),p);
% nodesOuter = find(abs(fdOuter([p]))<1e-3);
% nodesB = find(abs(fd([p]))<1e-3);
% 
% % Clear extra data
% clear('p');

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
