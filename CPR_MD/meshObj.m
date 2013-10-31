classdef meshObj
    % Mesh Class object
    % Create an unifor gird, set properties and elements.
    
    properties
        range
        numberElements
        solutionDegree
        elementSolPoints
    end
    
    methods
        function grid = meshObj(range,numberElements,PolyDeg)
            if nargin > 0 % Support calling with 0 arguments
                grid.range = range;
                grid.numberElements = numberElements;
                grid.solutionDegree = PolyDeg;
                grid.elementSolPoints = PolyDeg+1;
            end
        end
        function xxx
            
        end
    end
    
end

