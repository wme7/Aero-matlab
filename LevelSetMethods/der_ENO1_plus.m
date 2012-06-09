function [data_x] = der_ENO1_plus(data, dx)
%
% Calculates the derivative (plus) using
% first order accurate ENO scheme
% takes 1-D data
% data: input data
% dx: grid resolution
% Note: before entering this function, data needs to be 
% extended by 1 at the beginning and end (values don't matter)
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(1) = 2*data(2)-data(3);
data(end) = 2*data(end-1)-data(end-2);

data_x(2:end-1) = (data(3:end)-data(2:end-1))/dx;



