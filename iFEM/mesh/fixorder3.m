function [elem,idx,volume] = fixorder3(node,elem)
%% FIXORDER3 fix orientation of tetrahedron 
% 
%   elem = FIXORDER3(node,elem) computes signed volume of all tetrahedron
%   in the triangulation and switch the vertices such that all signed
%   volume is positive.
%   
%   [elem,idx,volume] = FIXORDER3(node,elem) also outputs the index set of
%   elements whose area is negative and the absolute value of area.
%
% See also fixorientation, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% compute signed volume of each tetrahedron
volume = simplexvolume(node,elem);
% find tetrahedron with negative volume and switch the vertices
idx = find(volume<0); 
elem(idx,[2 3]) = elem(idx,[3 2]);
volume(idx) = -volume(idx);