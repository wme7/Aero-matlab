
    %This example generates NE elements for a rectangular structure of
    %length = Ly units and width = Lx units with Nx divisions on the x   
    %axis using femTriangularMeshGenerator function
    
    %   Kehinde Orolu
    %   Systems Engineering
    %   University of Lagos, Nigeria
    %   olukeh@yahoo.com
    
    cla
    Lx=10;
    Ly=10;
    Nx=8;
    NE=144;
    
    [coords cT nNodes ]=femTriangularMeshGenerator(Lx,Ly,Nx,NE);
    
    disp(['Number of nodes =  ',num2str(nNodes)])
    disp('Connectivity Table')
    disp(cT)
    
    
    z=1;
    for i=1:NE
        figure(1),patch('Vertices',coords(z:z+2,:),'Faces',[1,2,3],'FaceColor','none','EdgeColor','g')
        hold on
        z=z+3;
    end
    
    figure(1),scatter(coords(:,1),coords(:,2),'MarkerFaceColor','r')
    
    hold off