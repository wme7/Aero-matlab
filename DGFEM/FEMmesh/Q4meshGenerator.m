function [VX VY EtoV nV nE]=Q4meshGenerator(Lx,Ly,xE,yE)
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
    fstats=@(p,t) fprintf('%d nodes, %d elements %.2f\n\n', ...
                      size(p,1),size(t,1));
    
    %% ELEMENT TO VERTICES CONNECTIVITY 
    % Build EtoV matrix
    nE = xE*yE;   % total number of elements
    nV = (xE+1)*(yE+1);  % total number of nodes
    eNodes = 4;
    
    EtoV = zeros(nE,eNodes); e=1;
    for i = 1:xE;
        for j = 1:yE;
            k = j+(yE+1)*(i-1);
            EtoV(e,:) = [k, (yE+1)+k, (yE+1)+k+1, k+1];
            e = e+1;
        end
    end

    %% COORDINATES GENERATION
    % build VX and VY arrays
    nx = linspace(Lx(1),Lx(2),xE+1);
    ny = linspace(Ly(1),Ly(2),yE+1);
    [x,y] = meshgrid(nx,ny);
    VX=zeros(nV,1); VY=zeros(nV,1); e=1;
    for i = 1:xE+1
        for j = 1:yE+1
            VX(e) = x(j,i);
            VY(e) = y(j,i);
            e = e+1; % element counter
        end
    end
    
    %% MESH STATUS
    fstats([VX,VY],EtoV);
end