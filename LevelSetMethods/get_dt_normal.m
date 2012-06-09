function [dt] = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs)
%
% Calculate the Euler time step.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

if alpha <= 0 | alpha >= 1 
    error('alpha needs to be between 0 and 1!');
end

maxs = max(max(H1_abs/dx + H2_abs/dy));
dt = alpha/(maxs+(maxs==0));

