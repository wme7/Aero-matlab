function [VX VY EtoV nV nE]=M2DmeshGenerator(mesh)
    %%  Triangular uniform mesh for a rectangular domain, using T3 shape 
    %   structure for finite element analysis
    %   by Manuel Diaz, NTU, 2013.11.28
    %
    %   INPUTS
    %   Lx:[a,b]=   width range of the rectangular structure
    %   Ly:[c,d]=   Height range of the rectangular structure
    %   xE      =   Number of elements on x- axis
    %   yE      =   Number of elements on y- axis
    %
    %   OUTPUTS
    %   VX      =   x coordinate each element vertex
    %   VY      =   y coordinate each element vertex
    %   EtoV	=   vertices connectivity
    %   nV      =   Number of Vertices
    %   nE      =   Number of elements

    % Define mesh status function
    fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));
    
    % Add mesh2d to search path
    addpath('E:\Documents\Matlab\aero-matlab\MeshGeneration\UnstructuredGrid\Mesh2d_v24');
              
    switch mesh
        case 1
            fprintf('Forward Step\n');
            % Grid (Polygon) Vertex
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

            % Grid Refiment function
            a = 0.6; b = 0.2; h = @(x,y) 1/35 + 1/5*sqrt( (x-a).^2+(y-b).^2 );

            % Element Size, hdata
            hdata = [];             % Initialize hdata class
            hdata.hmax = 1/2;   	% Domain element size.
            % hdata.edgeh=[2,1/20]; % Element size on specified geometry edges.
            hdata.fun = h;          % User defined size function.
            % Solver options
            % options.dhmax = 0.01;

            % Build mesh
            [vert,EtoV] = mesh2d(nodes,cnect,hdata);
        otherwise
            error('Case not available')
    end

    % MESH STATUS
    fstats(vert,EtoV);
    
    % Outputs
    VX=vert(:,1); VY=vert(:,2); nV=size(vert,1); nE=size(EtoV,1);
    
    %% Remove from search path
    rmpath('E:\Documents\Matlab\aero-matlab\MeshGeneration\UnstructuredGrid\Mesh2d_v24');
    
end