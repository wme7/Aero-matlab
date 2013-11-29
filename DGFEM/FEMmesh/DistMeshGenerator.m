function [VX VY EtoV nV nE]=DistMeshGenerator(mesh)
    %%  Triangular uniform mesh for a rectangular domain, using T3 shape 
    %   structure for finite element analysis
    %   by Manuel Diaz, NTU, 2013.11.28
    %
    %   INPUTS
    %   mesh    =   Select one of the precomputed cases
    %
    %   OUTPUTS
    %   VX      =   x coordinate each element vertex
    %   VY      =   y coordinate each element vertex
    %   EtoV	=   vertices connectivity
    %   nV      =   Number of Vertices
    %   nE      =   Number of elements

    % Add Distmesh to search path
    addpath('E:\Documents\Matlab\aero-matlab\MeshGeneration\UnstructuredGrid\distmesh')
    
    % Define mesh status function
    fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));
    
    switch mesh
        case 1
        fprintf('Rectangle with circular hole, refined at circle boundary\n');
        fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
        fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
        [p,EtoV]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
        fstats(p,EtoV);
       
        case 2
        fprintf('NACA0012 airfoil\n');
        hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4;
        a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];
        fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
        fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);
        fixx=1-htrail*cumsum(1.3.^(0:4)');
        fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
        fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
        box=[circx-circr,-circr; circx+circr,circr];
        h0=min([hlead,htrail,hmax]);
        [p,EtoV]=distmesh2d(fd,fh,h0,box,fix);
        fstats(p,EtoV);

        case 3
        fprintf('Forward Step\n');
        fd = inline('ddiff(drectangle(p,0,3,0,1),drectangle(p,0.6,3,0,0.2))','p');
        pfix = [0 0;0.6 0;0.6 0.2; 3 0.2;3 1;0 1];
        [p,EtoV] = distmesh2d(fd,@huniform,1/7,[0,0; 3,1],pfix);
        fstats(p,EtoV);
        
        otherwise
            error('case not defined');
    end

    % Prepare output
    VX=p(:,1); VY=p(:,2); nV=size(p,1); nE=size(EtoV,1);

    % Remove Distmesh from search path
    rmpath('E:\Documents\Matlab\aero-matlab\MeshGeneration\UnstructuredGrid\distmesh')