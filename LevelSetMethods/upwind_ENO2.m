function [data_x] = upwind_ENO2(data, F, dx)
%
% Calculates the upwind derivative using
% second order accurate ENO scheme
% takes 1-D data
% data: input data
% F: force field
% dx: grid resolution
% Note: before entering this function, data and F need to be 
% extended by 2 at the beginning and end (values don't matter)
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
data(2) = 2*data(3)-data(4);
data(1) = 2*data(2)-data(3);
data(end-1) = 2*data(end-2)-data(end-3);
data(end) = 2*data(end-1)-data(end-2);

%Generate the divided difference tables
%ignoring division by dx for efficiency
D1 = (data(2:end)-data(1:end-1));
% ignoring division by dx since this will cancel out
D2 = (D1(2:end)-D1(1:end-1))/2;
absD2 = abs(D2);

for i=1:(length(data)-4)
    if F(i+2) > 0  % use D-
        k = i-1;
    elseif F(i+2) < 0  % use D+
        k = i;
    else
        continue;
    end

    Q1p = D1(k+2); %D1k_half;

    if absD2(k+1) <= absD2(k+2) %|D2k| <= |D2kp1|
        c = D2(k+1); %D2k;
    else
        c = D2(k+2); %D2kp1;
    end

    % ignoring multiplication by dx since this will also cancel out
    Q2p = c*(2*(i-k)-1);

    data_x(i+2) = Q1p+Q2p;
    data_x(i+2) = F(i+2)*data_x(i+2)/dx;
end









