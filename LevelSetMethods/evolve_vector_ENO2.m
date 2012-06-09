function [delta] = evolve_vector_ENO2(phi, dx, dy, u_ext, v_ext)
%
% Finds the amount of evolution under a vector field
% based force and using 2nd order accurate ENO scheme
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

delta = zeros(size(phi)+4);
data_ext = zeros(size(phi)+4);
data_ext(3:end-2,3:end-2) = phi;
% first scan the rows
for i=1:size(phi,1)
    delta(i+2,:) = delta(i+2,:) + upwind_ENO2(data_ext(i+2,:), u_ext(i+2,:), dx);
end
% then scan the columns
for j=1:size(phi,2)
    delta(:,j+2) = delta(:,j+2) + upwind_ENO2(data_ext(:,j+2), v_ext(:,j+2), dy);
end
delta = delta(3:end-2,3:end-2);

