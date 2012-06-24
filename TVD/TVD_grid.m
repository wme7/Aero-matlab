function [x,y,Dimension] = TVD_grid(nx,ny)
% for nx and ny > 4, that means we would compute a 2D problem.
% for nx or  ny == 1, but not both, we would compute a 1D case.
% Compute the x and y vectors and Dimension of our problem

if nx >= 4 && ny >= 4
    x = 2:nx-1; 
    y = 2:ny-1;
    Dimension = 2;
elseif nx == 1 && ny >= 4
    x = 1;
    y = 2:ny-1;
    Dimension = 1;
elseif nx == 1 && ny <= 4
    error('Be careful!, nx == 1 && ny should be >= 4 for 1D problem')
elseif nx >= 4 && ny == 1
    x = 2:nx-1;
    y = 1;
    Dimension = 1;
elseif nx <= 4 && ny == 1
    error('Be careful!, nx should be >= 4 %% ny == 1 for 1D problem')
else
    error('Be careful!, nx and ny should not be to 1 at the same time!')
end

return