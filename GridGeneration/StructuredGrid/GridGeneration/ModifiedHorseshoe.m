function [x y] = ModifiedHorseshoe(Xi,Eta,R)

% R is aspect ratio of the major to minor axes of ellipse
% R >= 1

% Check the input arguments and set default values

if nargin == 2
    R = 3. ;
end

x = -(1+Eta)*cos(pi*Xi) ;
y = (1+(2*R-1)*Eta)*sin(pi*Xi) ;