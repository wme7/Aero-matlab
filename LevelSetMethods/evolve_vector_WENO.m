function [delta] = evolve_vector_WENO(phi, dx, dy, u_ext, v_ext)
%
% Finds the amount of evolution under a vector field
% based force and using 5th order accurate WENO scheme
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

delta = zeros(size(phi)+6);
data_ext = zeros(size(phi)+6);
data_ext(4:end-3,4:end-3) = phi;
% first scan the rows
for i=1:size(phi,1)
    delta(i+3,:) = delta(i+3,:) + upwind_WENO(data_ext(i+3,:), u_ext(i+3,:), dx);
end
% then scan the columns
for j=1:size(phi,2)
    delta(:,j+3) = delta(:,j+3) + upwind_WENO(data_ext(:,j+3), v_ext(:,j+3), dy);
end
delta = delta(4:end-3,4:end-3);

