function [data_x] = der_WENO_minus(data, dx)
%
% Calculates the derivative (minus) using
% fifth order accurate WENO scheme
% takes 1-D data
% data: input data
% dx: grid resolution
% Note: before entering this function, data needs to be 
% extended by 3 at the beginning and end (values don't matter)
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(3) = 2*data(4)-data(5);
data(2) = 2*data(3)-data(4);
data(1) = 2*data(2)-data(3);
data(end-2) = 2*data(end-3)-data(end-4);
data(end-1) = 2*data(end-2)-data(end-3);
data(end) = 2*data(end-1)-data(end-2);

D1 = (data(2:end)-data(1:end-1))/dx;

for i=1:(length(data)-6)
    k = i-1;
    v1 = D1(k+1);
    v2 = D1(k+2);
    v3 = D1(k+3);
    v4 = D1(k+4);
    v5 = D1(k+5);

        
    data_x_1 = v1/3 - 7*v2/6 + 11*v3/6;
    data_x_2 = -v2/6 + 5*v3/6 + v4/3;
    data_x_3 = v3/3 + 5*v4/6 - v5/6;

    epsilon = 1e-6*max([v1 v2 v3 v4 v5].^2) ...
        + 1e-99;

    S1 = (13/12)*(v1-2*v2+v3).^2 ...
        + 0.25*(v1-4*v2+3*v3).^2;
    S2 = (13/12)*(v2-2*v3+v4).^2 ...
        + 0.25*(v2-v4).^2;
    S3 = (13/12)*(v3-2*v4+v5).^2 ...
        + 0.25*(3*v3-4*v4+v5).^2;

    alpha1 = 0.1/((S1+epsilon).^2);
    alpha2 = 0.6/((S2+epsilon).^2);
    alpha3 = 0.3/((S3+epsilon).^2);
    a_total = alpha1+alpha2+alpha3;
    data_x(i+3) = (alpha1/a_total)*data_x_1 + (alpha2/a_total)*data_x_2 ...
        + (alpha3/a_total)*data_x_3;

    data_x(i+3) = data_x(i+3);
end

