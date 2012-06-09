function [data_x] = upwind_ENO1(data, F, dx)
%
% Calculates the upwind derivative using
% first order accurate ENO scheme
% takes 1-D data
% data: input data
% F: force field
% dx: grid resolution
% Note: before entering this function, data and F need to be 
% extended by 1 at the beginning and end (values don't matter)
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


if length(F) == 1
    F = F*ones(size(data));
end

if length(data) ~= length(F)
    error('force and data lengths does not match');
end

% if nargin < 3
%     alpha = 0.5;
% end

data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(1) = 2*data(2)-data(3);
data(end) = 2*data(end-1)-data(end-2);

%Generate the divided difference tables
D1 = (data(2:end)-data(1:end-1))/dx;

for i=1:(length(data)-2)
    if F(i+1) > 0  % use D-
        data_x(i+1) = F(i+1)*D1(i);
    elseif F(i+1) < 0  % use D+
        data_x(i+1) = F(i+1)*D1(i+1);
    else
        continue;
    end
end



