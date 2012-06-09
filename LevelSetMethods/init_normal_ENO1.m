function [Vn_ext] = init_normal_ENO1(Vn)
%
% Extends boundary of Vn
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

Vn_ext = zeros(size(Vn)+2);
Vn_ext(2:end-1,2:end-1) = Vn;

