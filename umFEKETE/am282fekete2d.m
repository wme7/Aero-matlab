function [X,Y] = am282fekete2d(umP) 

 % Copyright 2001, Brown University, Providence, Rhode Island.
 %
 % All Rights Reserved
 % 
 % Permission to use this software for noncommercial research and
 % educational purposes is hereby granted without fee.
 % Redistribution, sale, or incorporation of this software into a
 % commercial product is prohibited.
 % 
 % BROWN UNIVERSITY DISCLAIMS ANY AND ALL WARRANTIES WITH REGARD TO
 % THIS SOFTWARE,INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 % AND FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN
 % UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 % DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE,
 % DATA OR PROFITS.

 % Method taken from:
 %
 % "An Algorithm for computing Fekete Points in the triangle"
 % Mark Taylor, Beth Wingate, and R.E. Vincent. 
 % SIAM Journal on Numerical Analysis, Accepted 
 %
 % http://www.c3.lanl.gov/~wingate/papers.html
 % http://www.scd.ucar.edu/css/staff/taylor/points/fekete.html

 % Low storage Runge-Kutta coefficients
 rk4a = [  0.0   -567301805773.0/1357537059087.0  -2404267990393.0/2016746695238.0  -3550918686646.0/2091501179385.0  -1275806237668.0/842570457699.0];
 rk4b = [ 1432997174477.0/9575080441755.0   5161836677717.0/13612068292357.0  1720146321549.0/2090206949498.0  3134564353537.0/4481467310338.0  2277821191437.0/14882151754819.0];			     
 rk4c = [ 0.0  1432997174477.0/9575080441755.0  2526269341429.0/6820363962896.0  2006345519317.0/3224310063776.0  2802321613138.0/2924317926251.0 1.];

 % set the maximum polynomial order
 % umP = 6;

 % set the nodes to be equispaced
 umNpts = (umP+1)*(umP+2)/2;
 X = zeros(umNpts,1);
 Y = zeros(umNpts,1);

 skip = 1;
 for i=1:umP+1
  for j=1:umP+1 + (1-i)
   X(skip) = (-1*(umP-(i-1))+1*(i-1))/umP;
   Y(skip) = (-1*(umP-(j-1))+1*(j-1))/umP;
   skip = skip + 1;
  end
 end
 vmask  = [1,umP+1,umNpts];
 emask1 = [1:umP-1];
 emask2 = [1:umP-1];
 emask3 = [2:umP];
 skip = 1;
 for i=1:umP
  if(i>1)
   emask1(i-1) = skip;
  end
  for j=1:umP+1 + (1-i)
   skip = skip + 1;
  end
  if(i>1)
   emask2(i-1) = skip-1;
  end
 end

 Xorig = X;
 Yorig = Y;
 hold off;
 
 % compute initial Lebesgue constant
 am282lebesgue2d;

 % storage for node velocities
 U = zeros(umNpts,1);
 V = zeros(umNpts,1);

 % RK storage for X and Y time integrating residuals 
 resX = zeros(umNpts,1);
 resY = zeros(umNpts,1);

 % constants
 RKtime = 0;
 time   = 0;
 dt     = 0.5/(umP*umP);
 FinalTime = 1000;
 Nsteps = FinalTime/dt;

 % outer time step loop
 for tstep = 1:Nsteps

  % inner multi-stage Runge-Kutta loop
  for INTRK = 1:5
   % initiate Runge-Kutta stage
   resX   = rk4a(INTRK)*resX;
   resY   = rk4a(INTRK)*resY;
   RKtime = time+dt*rk4c(INTRK);

   % calculate derivative matrices
   [Dr,Ds] = am282derivmatrix2d(umP,X,Y);
 
   % set velocity of each node, to the gradient of 
   % its Lagrange interpolating polynomial 
   U = diag(Dr);
   V = diag(Ds);

   % constrain boundary nodes
   U(vmask)  = 0;
   V(vmask)  = 0;
   V(emask1) = 0;
   UdotN = U(emask2)+V(emask2);
   U(emask2) = U(emask2)-0.5*UdotN;
   V(emask2) = V(emask2)-0.5*UdotN;
   U(emask3) = 0;

   % update R-K residuals
   resX = resX+dt*U;
   resY = resY+dt*V;

   % finish Runge-Kutta stage
   X = X+rk4b(INTRK)*resX;
   Y = Y+rk4b(INTRK)*resY;
  end;
  % end inner-RK loop

  if(tstep == 0 | mod(tstep,1) == 0)
   % plot original node positions

   hold off;
   % calculate/plot Lebesgue function, and Lebesgue number
   am282lebesgue2d;
   
   hold on;
   % plot original nodes
   scatter3(Xorig,Yorig,-5*ones(size(X)),'r');
   % plot current node positions
   scatter3(X,Y,-5*ones(size(X)),'g');
   axis manual; axis([-1,1,-1,1,-6,1.5*umP]);

   hold off; pause(.05); drawnow;
 
  end

  % update time
  time = time+dt;

  % stopping criterion
  maxU = max(abs(U));
  maxV = max(abs(V));
  maxVelocity = [maxU,maxV]
  if( (maxU < 1e-10/dt) & (maxV < 1e-10/dt))
   break;
  end
 end
