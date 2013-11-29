function [VX VY EtoV nV nE]=AGmeshGenerator(Lx,Ly,xE,yE,coordinate)
    %%  Triangular uniform mesh for a rectangular domain, using T3 shape 
    %   structure for finite element analysis
    %   by Manuel Diaz, NTU, 2013.11.28
    %
    %   INPUTS
    %   Lx:[a,b]=   width range of the rectangular structure
    %   Ly:[c,d]=   Height range of the rectangular structure
    %   xE      =   Number of elements on x- axis
    %   yE      =   Number of elements on y- axis
    %   Coordinate = One of the available coordinate system
    %
    %   OUTPUTS
    %   VX      =   x coordinate each element vertex
    %   VY      =   y coordinate each element vertex
    %   EtoV	=   vertices connectivity
    %   nV      =   Number of Vertices
    %   nE      =   Number of elements
    
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

    % Coordinates systems available:
    %   PolarCoordinates
    %   ParabolicCylinderCoordinates
    %   EllipticCylinderCoordinates
    %   Horseshoe
    %   ModifiedHorseshoe
    %   BipolarCoordinates
    
    % discretize along xi and eta axis
    xi = linspace(Lx(1),Lx(2),xE+1);
    eta = linspace(Ly(1),Ly(2),yE+1);
    [XI,ETA] = meshgrid(xi,eta);
    vert = AnalyticGrids(XI,ETA,coordinate);
    X = vert.x; Y = vert.y;
    % Build VX and VY arrays
    VX=zeros(nV,1); VY=zeros(nV,1); e=1;
    for i = 1:xE+1
        for j = 1:yE+1
            VX(e) = X(j,i);
            VY(e) = Y(j,i);
            e = e+1; % element counter
        end
    end  
       
end