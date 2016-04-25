function [idxs, delta] = evalDelta( Xl, X, Y, N, h);
%
% [idxs delta] = evalDelta( Xl, X, Y, N, h );
%
%  Calculates the delta function at grid locations, idxs
%
%  Returns:
%     idxs  = L by 16 matrix, where L is the number of Langrangian Pts
%             gives the indices in the Euclidean mesh of where
%             the weight should be applied
%     delta = the weights to apply at each index above
%     
%  Input:
%     Xl    = L by 2 array of the location of the Langrangian Pts
%     X     = N*N length vector of all X spatial locations
%     Y     = " but in Y direction
%     N     = number of mesh points in each direction
%     h     = mesh width
%
%  NOTES:
%     MESH POINTS ARE ASSUME TO BE AT LOCATIONS (i,j)*h, i=0..N-1,
%     j=0..N-1. Does not handle wrapping across periodic boundaries.
%
%
%
%  License: This code is free to use for any purposes, provided
%           any publications resulting from the use of this code
%           reference the original code/author.
%
%  Author:  Samuel Isaacson (isaacson@math.utah.edu)
%  Date:    11/2007
%
%  Please notify the author of any bugs, and contribute any
%  modifications or bug fixes back to the original author.
%
%  Disclaimer:
%   This code is provided as is. The author takes no responsibility 
%   for its results or effects.


persistent isInit xL xR yD yU idxVecX idxVecY;

if( size( isInit ) == 0 )
  isInit = 1;
  
  xL = -2:1;
  xR = -1:2;
  yD = [-2:1];
  yU = [-1:2];
  
  idxVecX = zeros(1, 16);
  idxVecY = zeros(1, 16);
end


XlLen = length(Xl(:,1));
idxs  = zeros(XlLen, 16);
delta = zeros(XlLen, 16);

for( l = 1:XlLen )
  
  % cell Xl lies in, assume Xmesh = j * h:
  fIdx    = Xl(l,:) / h;
  cellIdx = round( fIdx );

  % find quadrant of this cell, and hence mesh indices to use, note
  % indices start at zero not one!:
  if( cellIdx(1) > fIdx(1) )
    xIdx = xL + cellIdx(1);
  else
    xIdx = xR + cellIdx(1);    
  end
  
  if( cellIdx(2) > fIdx(2) )
    yIdx = yD + cellIdx(2);    
  else
    yIdx = yU + cellIdx(2);
  end
  
  idxVecX = [ xIdx xIdx xIdx xIdx ]';
  idxVecY = [yIdx(1); yIdx(1); yIdx(1); yIdx(1); yIdx(2); yIdx(2); yIdx(2); yIdx(2); yIdx(3);  yIdx(3);  yIdx(3); yIdx(3); yIdx(4); yIdx(4); yIdx(4); yIdx(4);];
  
  
  delta(l,:) = (evalPhi( fIdx(1) - idxVecX ) .* evalPhi( fIdx(2) - idxVecY ) )' /(h*h);
  
  xIdx = mod(xIdx, N) + 1;
  yIdx = mod(yIdx, N);

  idxs(l,:) = [ (xIdx + N * yIdx(1)) (xIdx + N * yIdx(2)) (xIdx + N * yIdx(3)) (xIdx + N * yIdx(4))];  
  
end