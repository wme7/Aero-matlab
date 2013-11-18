function [x y] = PolarCoordinates(Xi,Eta,theta0,theta1,r0,r1) 
% theta0, theta1 are two angles. 
% 0<=theta0<theta1<=2*pi

% r0, r1 are two radii
% 0<=r0<r1

% Check for the input arguments and set default values

if nargin == 2
    theta0 = 0. ;
    theta1 = -pi/4 ;

    r0 = 1. ;
    r1 = 2. ;
    
end

theta = (theta0-theta1)*Xi ;
r = r0+(r1-r0)*Eta ;

x = r*cos(theta) ;
y = r*sin(theta) ;