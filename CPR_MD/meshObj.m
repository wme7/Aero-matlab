classdef meshObj
    % Mesh Class object
    % Create an unifor gird, set properties and elements.
    
    properties
        range
        numberElements
        solutionDegree
        elementSolPoints
        mesh
    end
    
    methods
        function obj = meshObj(range,numberElements,PolyDeg)
            if nargin > 0 % Support calling with 0 arguments
                obj.range = range;
                obj.numberElements = numberElements;
                obj.solutionDegree = PolyDeg;
                obj.elementSolPoints = PolyDeg+1;
            end
        end
        function uniformMesh(obj)
            obj.mesh = obj.range(1):0.1:obj.range(2);
        end
    end
    
end

