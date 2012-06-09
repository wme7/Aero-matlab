function [dt] = get_dt_vector(alpha,dx,dy,u,v)
%
% Calculate the Euler time step.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

if alpha <= 0 | alpha >= 1 
    error('alpha needs to be between 0 and 1!');
end

maxs = max(abs(u(:))/dx + abs(v(:))/dy);
dt = alpha/(maxs+(maxs==0));





