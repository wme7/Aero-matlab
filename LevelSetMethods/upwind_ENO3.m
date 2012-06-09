function [data_x] = upwind_ENO3(data, F, dx)
%
% Calculates the upwind derivative using
% third order accurate ENO scheme
% takes 1-D data
% data: input data
% F: force field
% dx: grid resolution
% Note: before entering this function, data and F need to be 
% extended by 3 at the beginning and end (values don't matter)
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
data(3) = 2*data(4)-data(5);
data(2) = 2*data(3)-data(4);
data(1) = 2*data(2)-data(3);
data(end-2) = 2*data(end-3)-data(end-4);
data(end-1) = 2*data(end-2)-data(end-3);
data(end) = 2*data(end-1)-data(end-2);

%Generate the divided difference tables
%ignoring division by dx for efficiency
D1 = (data(2:end)-data(1:end-1));
D2 = (D1(2:end)-D1(1:end-1))/2;
absD2 = abs(D2);
D3 = (D2(2:end)-D2(1:end-1))/3;
absD3 = abs(D3);

for i=1:(length(data)-6)
    if F(i+3) > 0  % use D-
        k = i-1;
    elseif F(i+3) < 0  % use D+
        k = i;
    else
        continue;
    end

    Q1p = D1(k+3); %D1k_half;

    if absD2(k+2) <= absD2(k+3) %|D2k| <= |D2kp1|
        kstar = k-1;
        c = D2(k+2);
    else
        kstar = k;
        c = D2(k+3);
    end
    Q2p = c*(2*(i-k)-1);

    if absD3(kstar+2) <= absD3(kstar+3) %|D3kstar_half| <= |D3kstar_1_half|
        cstar = D3(kstar+2); %D3kstar_half;
    else
        cstar = D3(kstar+3); %D3kstar_1_half;
    end
    Q3p = cstar*( 3*(i-kstar)*(i-kstar) - 6*(i-kstar) + 2 );

    data_x(i+3) = Q1p+Q2p+Q3p;

    data_x(i+3) = F(i+3)*data_x(i+3)/dx;
end





