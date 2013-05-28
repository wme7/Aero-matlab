function [nodeCoordinates elementNodes] = Mesh_Lab10(name)
% Mesh for Lab10 problems

switch name
    case 'T3' % Using T3 elements
        
        nodeCoordinates = [0 0;0 10;10 0;10 10;20 0;20 10;30 0;30 10];
        elementNodes =  [1 3 4;1 4 2;3 5 6;3 6 4;5 7 8;5 8 6];
        drawingMesh(nodeCoordinates,elementNodes,'T3','b-o');
        
    case 'Q4' % Using Q4 elements
        
        nodeCoordinates = [0 0;0 10;10 0;10 10;20 0;20 10;30 0;30 10];
        elementNodes =  [1 3 4 2;3 5 6 4;5 7 8 6];
        drawingMesh(nodeCoordinates,elementNodes,'Q4','b-o');
        
    case 'Q8' % Using Q8 elements
        
        nodeCoordinates = [0 0;0 5;0 10; ...
                            5 0; 5 10; ...
                           10 0;10 5;10 10; ...
                            15 0;15 10; ...
                           20 0;20 5;20 10; ...
                            25 0;25 10; ...
                           30 0;30 5;30 10];
        elementNodes =  [1 6 8 3 4 7 5 2; ...
                         6 11 13 8 9 12 10 7; ...
                         11 16 18 13 14 17 15 12];
        %drawingMesh(nodeCoordinates,elementNodes,'Q8','b-o');
        scatter(nodeCoordinates(:,1),nodeCoordinates(:,2))
        
    case 'Q12' % Using Q12 elements
        
        nodeCoordinates = [0 0;0 10/3;0 20/3;0 10;
                            10/3 0;10/3 10;
                            20/3 0;20/3 10;
                           10 0;10 10/3;10 20/3;10 10;
                            10+10/3 0;10+10/3 10;
                            10+20/3 0;10+20/3 10;
                           20 0;20 10/3;20 20/3;20 10;
                            20+10/3 0;20+10/3 10;
                            20+20/3 0;20+20/3 10;
                           30 0;30 10/3;30 20/3;30 10];
        elementNodes =  [1 9 12 4 5 7 10 11 8 6 3 2;
                         9 17 20 12 13 15 18 19 16 14 11 10;
                         17 25 28 20 21 23 26 27 24 22 19 18];
        %drawingMesh(nodeCoordinates,elementNodes,'Q12','b-o');
        scatter(nodeCoordinates(:,1),nodeCoordinates(:,2))
        
    otherwise
        error('wrong case');
end