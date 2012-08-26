  % Math578 - HW4 - Discontinuous Galerkin (with Legendre pols)
  % for the advection equation q_t + u*q_x=0

  % DESCRIBE THE ORDER OF THE SCHEME
  % Order of polynomials used for interpolation
  p = 4;

  % Chebychev nodes on the interval [-1,1]
  xcheby = sin(0.5*pi*linspace(-1,1,p+1)');

  % DIVIDE THE PIPE INTO CELLS
  % number of cells to divide interval into
  Ncells = 100;

  % Physical limits of the pipe
  xmin = -4;
  xmax =  4;

  % set Number of nodes along pipe
  Nnodes = Ncells+1;
  nodex = linspace(xmin,xmax,Nnodes);

  % for each cell of the pipe fill with p+1 Chebychev nodes
  x = zeros(p+1, Ncells);
  for i=1:Ncells
    x(:,i) = ((nodex(i+1)-nodex(i))/2)*xcheby + (nodex(i)+nodex(i+1))/2;
  end

  % CREATE THE MATRICES NEEDED FOR THE DG SCHEME
  % Create the Vandermonde matrix for the Chebychev nodes
  V = LEGvdm(xcheby,p);

  %Differentiation matrix 
  %(acts on the vector of coefficients of q 
  % and gives coefficients of q_x)

  Dhat = zeros(p+1);
  
  for n=0:p
    for m=n+1:2:p
      Dhat(n+1,m+1) =  (2*n+1);
    end
  end

  M = diag(2./(2*(0:p)+1));

  % Stiffness matrix
  D = M*Dhat;

  % Modify D
  % i.e. D <- D(n,m)+ 2*(-1)^(n+m); 0<=m,n<=p

  Dmod = zeros(p+1);
  for i = 1:p+1
    for j = 1:p+1
      Dmod(i,j) = D(i,j) + (-1)^((i-1)+(j-1));
    end
  end


  % surface term matrix
  S = (-1).^[0:p]'* ones(1,p+1);
	
  % INITIAL CONDITIONS
  %Initial condition q in each cell and its interpol. coeffs ro
  %(cells indexed by i, each cell has p+1 nodes)
  q   = zeros(p+1, Ncells);
  
  q(:)   = exp(-x(:).^2);
  rho = V\q;

  % PARAMETERS FOR THE SCHEME
  u = 1.;

  % width of each cell 
  dx = nodex(2:Nnodes)-nodex(1:Nnodes-1);

  % time step for time marching
  dt = min(min(dx))/((p+1)^2);

  % final time to integrate advection equation to
  T = 4;                  

  % number of time steps
  Ntsteps = floor(T/dt); 

  % coefficients 
  coeff   = (2*(0:p)+1)'*(u./(dx));

  % TIME STEPPING

  for tstep = 1:Ntsteps

    sigma = rho;

    for rkstage = 3:-1:1
      % stage 1
      volterms = -Dmod*sigma;
      
      % flux coming in from left neighbors
      surfterm =  S*sigma;
      
      idthis = 2:Ncells;
      idneig = 1:Ncells-1;
    
      volterms(:, idthis) = ...
	  (dt/rkstage)*coeff(:,idthis).*(volterms(:,idthis)+surfterm(:,idneig));
    
      volterms(:, 1) = (dt/rkstage)*coeff(:,1).*volterms(:,1);

      sigma    = rho + volterms;
    end
    
    rho = sigma;

    if ( ~mod(tstep, 100)  )
      q = V*rho;

%      plot(x(:),abs(q(:)-exp(-(x(:)-dt*tstep).^2)));
      plot(x(:),q(:));
      pause(0.1);
    end
  end
