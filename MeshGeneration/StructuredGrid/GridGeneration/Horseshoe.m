function [x y] = Horseshoe(Xi,Eta,row,b0,b1)

% row is aspect ratio of the major to minor axes of ellipse
% row > 0

% b0, b1 are two radii
% 0<b0<b1

% check for the input arguments and set default values
if nargin == 2
    row = 3. ;
end

if nargin == 2
    b0 = 4. ;
    b1 = 7. ;
end

r = b0+(b1-b0)*Eta ;
theta = pi/2*(1-2*Xi) ;

x = row*r*cos(theta) ;
y = r*sin(theta) ;
