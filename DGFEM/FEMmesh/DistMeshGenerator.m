function [VX VY EtoV nV nE]=DistMeshGenerator()
    %%  Triangular uniform mesh for a rectangular domain, using T3 shape 
    %   structure for finite element analysis
    %   by Manuel Diaz, NTU, 2013.11.28
    %
    %   OUTPUTS
    %   VX      =   x coordinate each element vertex
    %   VY      =   y coordinate each element vertex
    %   EtoV	=   vertices connectivity
    %   nV      =   Number of Vertices
    %   nE      =   Number of elements

fd = inline('ddiff(drectangle(p,0,3,0,1),drectangle(p,0.6,3,0,0.2))','p');
pfix = [0 0;0.6 0;0.6 0.2; 3 0.2;3 1;0 1];
[vert,EtoV] = distmesh2d(fd,@huniform,1/7,[0,0; 3,1],pfix);
VX=vert(:,1); VY=vert(:,2); nV=size(vert,1); nE=size(EtoV,1);