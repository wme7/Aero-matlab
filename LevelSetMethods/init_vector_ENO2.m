function [u_ext, v_ext] = init_vector_ENO2(u,v)
%
% Extends boundaries of u and v
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

u_ext = zeros(size(u)+4);
v_ext = zeros(size(v)+4);
u_ext(3:end-2,3:end-2) = u;
v_ext(3:end-2,3:end-2) = v;


