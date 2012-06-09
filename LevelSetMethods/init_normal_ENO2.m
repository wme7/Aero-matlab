function [Vn_ext] = init_normal_ENO2(Vn)
%
% Extends boundary of Vn
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

Vn_ext = zeros(size(Vn)+4);
Vn_ext(3:end-2,3:end-2) = Vn;

