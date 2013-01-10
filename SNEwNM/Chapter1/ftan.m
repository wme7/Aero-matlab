function [y,yp] = ftan(x)
% FTAN   Example function with optional derivative.
%        [Y,YP] = FATAN(X) returns Y=tan(X) - X and
%        (optionally) YP = sec(X)^2 - 1
%
y = tan(x)-x;
if nargout == 2
    yp = sec(x)^2-1;
end
