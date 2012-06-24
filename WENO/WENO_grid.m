function [x,y,Dimension] = WENO_grid(nx,ny)
% for nx and ny > 7, that means we would compute a 2D problem.
% for nx or  ny == 1, but not both, we would compute a 1D case.
% Compute the x and y vectors and Dimension of our problem

if nx >= 7 && ny >= 7
    x = 4:nx-3; 
    y = 4:ny-3;
    Dimension = 2;
elseif nx == 1 && ny >= 7
    x = 1;
    y = 4:ny-3;
    Dimension = 1;
elseif nx == 1 && ny <= 7
    error('Be careful!, nx == 1 && ny should be >= 7 for 1D problem')
elseif nx >= 7 && ny == 1
    x = 4:nx-3;
    y = 1;
    Dimension = 1;
elseif nx <= 7  && ny == 1
    error('Be careful!, nx should be >= 7 %% ny == 1 for 1D problem')
else
    error('Be careful!, nx and ny should not be to 1 at the same time!')
end

return