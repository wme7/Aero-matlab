function [u_ext, v_ext] = init_vector_ENO1(u,v)
%
% Extends boundaries of u and v
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

u_ext = zeros(size(u)+2);
v_ext = zeros(size(v)+2);
u_ext(2:end-1,2:end-1) = u;
v_ext(2:end-1,2:end-1) = v;


