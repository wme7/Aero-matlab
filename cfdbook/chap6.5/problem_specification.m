%PROBLEM SPECIFICATION
% Specification of parameters, grid and boundary conditions

global Re

tic
if geval == 1	% Horizontal Poiseuille flow to the right
  dt = 0.1;	% Time step
  tend = 10;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = 0.75;	% Estimates of max. velocity components required for stability
  vmax = 0;	% condition for Adams-Bashforth-crank-Nicolson scheme
  
  xseglen = [1,1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [4,4];	% Number of cells along x-segments.
  ny = [4,4];	% Number of cells along y-segments.
		
  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [1,1; 1,1];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [2,2; 3,3];
    
elseif geval == 2	% Vertical Poiseuille flow upward
  dt = 0.5;	% Time step
  tend = 10;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = 0;	% Estimates of max. velocity components required for stability
  vmax = 0.75;	% condition for Adams-Bashforth-crank-Nicolson scheme

  xseglen = [1,1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [8,16];	% Number of cells along x-segments.
  ny = [8,8];	% Number of cells along y-segments.
		
  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [2,2; 3,3];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [1,1; 1,1];
  
elseif geval == 3	% Backward facing step
  dt = 3;	% Time step
  tend = 140;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = 1.25;	% Estimates of max. velocity components required for stability
  vmax = 0.25;	% condition for Adams-Bashforth-Crank-Nicolson scheme
  		% (not used in program ns1)

  xseglen = [5,5];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [20,10];	% Number of cells along x-segments.
  ny = [25,25];	% Number of cells along y-segments.

  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [1,1; 1,1];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [1,2; 3,3];
  
elseif geval == 4	% Driven cavity
  dt = 0.5;	% Time step
  tend = 25;	% End time
  Re = 200;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 0.55;  % Omega-scheme for time-stepping
  umax = 1;	% Estimates of max. velocity components required for stability
  vmax = 0;	% condition for Adams-Bashforth-crank-Nicolson scheme

  xseglen = [1,1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [15,15];	% Number of cells along x-segments.
  ny = [15,15];	% Number of cells along y-segments.
  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [1,1; 2,2];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [1,1; 1,1];
  
elseif geval == 5	% Uniform flow under angle alpha in [0,pi/2] without
			% segmentation
  alpha = pi/4;
  dt = 3;	% Time step
  tend = 60;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = abs(cos(alpha));% Estimates of max. velocity components required for stability
  vmax = abs(sin(alpha));% condition for Adams-Bashforth-crank-Nicolson scheme

  xseglen = [1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1];	% Similar to xseglen in y-direction.
  nx = [2];	% Number of cells along x-segments.
  ny = [2];	% Number of cells along y-segments.

  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [2; 3];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [2; 3];
  
elseif geval == 6	% Uniform flow under angle alpha in [0,pi/2] with
			% segmentation
  alpha = pi/4;
  dt = 3;	% Time step
  tend = 60;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = abs(cos(alpha));% Estimates of max. velocity components required for stability
  vmax = abs(sin(alpha));% condition for Adams-Bashforth-crank-Nicolson scheme

  xseglen = [1,1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [2,6];	% Number of cells along x-segments.
  ny = [2,3];	% Number of cells along y-segments.

  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [2,2; 3,3];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [2,2; 3,3];
  
elseif geval == 7	% Horizontal Poiseuille flow to the left
  dt = 0.5;	% Time step
  tend = 60;	% End time
  Re = 100;	% Reynolds number based on unit length and unit velocity
  central = 1;	% Enter 1 for central scheme or something else for upwind scheme
  omega = 1/2;  % Omega-scheme for time-stepping
  umax = 0.75;	% Estimates of max. velocity components required for stability
  vmax = 0;	% condition for Adams-Bashforth-crank-Nicolson scheme

  xseglen = [1,1];	% Enter segment lengths of horizontal boundary in
		% order of increasing x. Number of segments is arbitrary.
		% These segments are used for grid generation and boundary 
		% conditions. 
  yseglen = [1,1];	% Similar to xseglen in y-direction.
  nx = [4,4];	% Number of cells along x-segments.
  ny = [4,4];	% Number of cells along y-segments.
		
  % Type of boundary condition: 1: no-slip
  %				2: inflow
  %				3: outflow		
  % Give type of boundary condition along lower (first row) and upper 
  %(second row) horizontal segments in order of increasing x:
   
  xbc = [1,1; 1,1];
   	
  % Give type of boundary condition along left (first row) and right 
  % (second row) vertical segments in order of increasing y:
  
  ybc = [3,3; 2,2];    
  
else
  error('Wrong value in input for parameter geval')  	
end
		
if (size(xseglen)~=size(nx))|(size(yseglen)~=size(ny))
  error('Wrong correspondence between xseglen,nx or yseglen,ny in problem_specification')
end
if (size(xbc(1,:))~=size(nx))|(size(ybc(1,:))~=size(nx))
  error('Wrong correspondence between nx,xbc or ny,ybc in problem_specification')
end

tijd = toc; disp(['problem_specification time = ',num2str(tijd)])
