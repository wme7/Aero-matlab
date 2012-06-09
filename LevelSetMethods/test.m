%
% This is a test m-file to try out this package on 
% level set methods. It creates a signed distance function
% whose zero level set (the curve) is a square. In the default setup,
% this square is then shrunk under a force in the normal direction.
% Try using a combination of normal, vector field-based and 
% curvature-based forces to evolve the curve. 
% evolve2D() is the core function that handles the evolution.
% Type help evolve2D to see the documentation of this function.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%


phi = ones(100,100);
phi(30:60,30:60) = -1;
phi = double((phi > 0).*(bwdist(phi < 0)-0.5) - (phi < 0).*(bwdist(phi > 0)-0.5));

% Using reinit_SD on the step function instead may 
% smooth out the curve corners slightly (not desired).
% phi = reinit_SD(phi, dx, dy, 0.5, 'WENO', 200);

% External vector field
u = -0.5*(1/sqrt(2)).*ones(100, 100); %x component of the vector field
v = -0.5*(1/sqrt(2)).*ones(100, 100); %y component of the vector field

% b: weighting for curvature-based force. b needs to be positive
b = 0.3*ones(100, 100);

% Vn: force in the normal direction
Vn = -0.2*ones(100, 100);

dx=1;
dy=1;

figure;contour(phi,[0 0],'b');
phi = evolve2D(phi,dx,dy,0.5,15,'ENO3',1,1,Vn,0,u,v,0,b);
hold on; contour(phi,[0 0],'r'); hold off; axis equal;


