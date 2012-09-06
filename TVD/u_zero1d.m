function u_0 = u_zero1d(x)
%% Initial Condition (IC)
% This subroutine creates an special IC for testing purposes.
% The present IC models:
%
% 1. A square jump condition, from 0 to 1
% 2. A triangular jump condition, from 0 to 1
% 3. A sinusoidal jump condition, from 0 to 1
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06

%% IC
% number of points required
 n = length(x);

% Preallocating
u_0 = zeros(1,n); 

% Reference points
x0 = 1;             % point 0
x1 = ceil(1*n/10);  % point 1
x2 = ceil(3*n/10);  % point 2
x3 = ceil(4*n/10);  % point 3
x3_5 = ceil(5*n/10);% between 3 and 4 
x4 = ceil(6*n/10);  % point 4
x5 = ceil(7*n/10);  % point 5
x6 = ceil(9*n/10);  % point 6
x7 = n;             % point 7

% y references
y1 = 1;
y2 = 2;

% Load regions 
for x = 1:n
    if x0<=x && x<x1
    u_0(x) = y1;     % reference region
    elseif x1<=x && x<x2
    u_0(x) = y2;                            % ones region
    elseif x2<=x && x<x3
    u_0(x) = y1;     % reference region
    elseif x3<=x && x<x3_5
    u_0(x) = (y2-y1)/(x3_5-x3)*(x-x3)+y1;   % half triangle
    elseif x3_5<=x && x<x4
    u_0(x) = (y1-y2)/(x4-x3_5)*(x-x3_5)+y2; % half triangle
    elseif x4<=x && x<x5
	u_0(x) = y1;     % reference region
    elseif x5<=x && x<x6
	u_0(x) = sin((x-x5)/(x6-x5)*pi)+y1;      % semicircle region
    elseif x6<=x && x<=x7
	u_0(x) = y1;     % reference region
    else
    u_0(x) = y1;     % reference region
    end
end
