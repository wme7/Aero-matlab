function [delta] = evolve_vector_ENO1(phi, dx, dy, u_ext, v_ext)
%
% Finds the amount of evolution under a vector field
% based force and using 1st order accurate ENO scheme
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


delta = zeros(size(phi)+2);
data_ext = zeros(size(phi)+2);
data_ext(2:end-1,2:end-1) = phi;
% first scan the rows
for i=1:size(phi,1)
    delta(i+1,:) = delta(i+1,:) + upwind_ENO1(data_ext(i+1,:), u_ext(i+1,:), dx);
end
% then scan the columns
for j=1:size(phi,2)
    delta(:,j+1) = delta(:,j+1) + upwind_ENO1(data_ext(:,j+1), v_ext(:,j+1), dy);
end
delta = delta(2:end-1,2:end-1);

