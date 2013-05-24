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
        
        nodeCoordinates = [0 0;0 10;10 0;10 10;20 0;20 10;30 0;30 10];
        elementNodes =  [1 3 4 2;3 5 6 4;5 7 8 6];
        drawingMesh(nodeCoordinates,elementNodes,'Q4','b-o');
        
    case 'Q12' % Using Q12 elements
        
        nodeCoordinates = [0 0;0 10;10 0;10 10;20 0;20 10;30 0;30 10];
        elementNodes =  [1 3 4 2;3 5 6 4;5 7 8 6];
        drawingMesh(nodeCoordinates,elementNodes,'Q4','b-o');
        
    otherwise error('wrong case');
        
end