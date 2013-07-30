function [coords cT nNodes ]=femTriangularMeshGenerator(Lx,Ly,Nx,NE)
    %   This function generates triangular mesh for a rectangular
    %   shape structure for finite element analysis
    %   [coords cT nNodes ]=femTriangularMeshGenerator(Lx,Ly,Nx,NE)
    %   coords  =   x and y coordinates of each element nodes
    %   cT      =   nodal connectivity
    %   nNodes  =   Number of nodes
    %   Lx      =   width of the rectangular structure
    %   Ly      =   Height of the rectangular structure
    %   Nx      =   Number of divisions on x- axis
    %   NE      =   Number of elements
    %   
    
    %   Kehinde Orolu
    %   Systems Engineering
    %   University of Lagos, Nigeria
    %   olukeh@yahoo.com
    
    if mod((NE/Nx),2)~=0
        errordlg('The No of divisions on X axis must divide No of Elements twice')
        
    end
    
    Ny=NE/(2*Nx);   %Divisions on y axis
    
    nNodes =(Nx+1)*(Ny+1);  %No of nodes
    
    m=1;
    j=1:Nx;
    k=linspace(Nx*2,NE,Ny);
    
    
    
    for i=1:Ny
        cT(m:2:k(i),1)= j;  %node 1 of 1st element
        cT(m+1:2:k(i),1)= j;  %node 1 of 2nd element
        cT(m:2:k(i),2)=j+1; %%node 2 of 1st element
        cT(m+1:2:k(i),2)=j+Nx+2;%node 2 of 2nd element
        cT(m:2:k(i),3)= j+Nx+2;  %node 3 of 1st element
        cT(m+1:2:k(i),3)=j+1+Nx;    %%node 3 of 1st element
        
        m=k(i)+1;
        j=j+Nx+1;
    end
    %%%%%%%%%%%%%%%COORDINATES GENERATION%%%%%%%%%%%%%%%%%%%
    
    ax=linspace(0,Lx,Nx+1); %%%x coordinates
    by=linspace(0,Ly,Ny+1); %%%y coordinates
    X1=[];
    Y1=[];
    for i1=1:Ny+1
        %    General Nodal Coordinates layer by layer
        by1(1:Nx+1)=by(i1);
        X1=[X1 ax];
        Y1=[Y1 by1];
    end
    j=1:3;
    
    %each element coordinates
    for n=1:NE
        X(j,1) = X1(cT(n,:));
        Y(j,1)=Y1(cT(n,:));
        j=j+3;
    end
    coords=[X Y];   %x and y coordinates