function setgrid(grid,range,numberElements,PolyDeg)
    if nargin > 0 % Support calling with 0 arguments
        grid.range = range;
        grid.numberElements = numberElements;
        grid.SolutionDegree = PolyDeg;
        grid.elementSolPoints = PolyDeg+1;
    end
end