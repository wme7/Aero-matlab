function [y,yp] = fatan(x)
% FATAN  Arctangent function with optional derivative.
%        [Y,YP] = FATAN(X) returns Y=atan(X) and
%        (optionally) YP = 1/(1+X^2)
%
%
y = atan(x);
if nargout == 2
    yp = 1/(1+x^2);
end
