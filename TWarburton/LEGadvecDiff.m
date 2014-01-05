% Math578 - HW4 - Discontinuous Galerkin (with Legendre pols)
% for the advection diffusion equation q_t + u*q_x= Dcoeff*q_xx

% DESCRIBE THE ORDER OF THE SCHEME
% Order of polynomials used n interpolation
p = 4;

% Chebychev nodes on the interval [-1,1]
xcheby = sin(0.5*pi*linspace(-1,1,p+1)');

% DIVIDE THE PIPE INTO CELLS
% number of cells to divide interval into
Ncells = 20;

% Physical limits of the pipe
xmin = 0;
xmax = 2*pi;

% set Number of nodes along pipe
Nnodes = Ncells+1;
nodex = linspace(xmin,xmax,Nnodes);

% in each cell of the pipe fill with p+1 Chebychev nodes
x = zeros(p+1, Ncells);
for i=1:Ncells
  x(:,i) = ((nodex(i+1)-nodex(i))/2)*xcheby + (nodex(i)+nodex(i+1))/2;
end

% CREATE THE MATRICES NEEDED IN THE DG SCHEME
% Create the Vandermonde matrix with Chebychev nodes
V = LEGvdm(xcheby,p);

%Differentiation matrix 
Dhat = zeros(p+1);
for n=0:p
  for m=n+1:2:p
    Dhat(n+1,m+1) =  (2*n+1);
  end
end

% M = int_{-1}^{1} L_n(x)*L_m(x) dx (note dx=1 assumed)
M = diag(2./(2*(0:p)+1));

% D_nm = int_{-1}^{1}  L_n(x)*L_m(x)*dx
D = M*Dhat;

% F_nm = L_n(-1)*L_m(-1)
F = (-1).^[0:p]'*(-1).^[0:p];

% G_nm = -L_n(-1)*L_m(1)
G = -(-1).^[0:p]'*ones(1,p+1);

% H_nm =  -L_n(1)*L_m(1)
H = -ones(p+1,p+1);

% J_nm =  L_n(1)*L_m(-1)
J = ones(p+1,1)*(-1).^[0:p];

% multiply all the matrices by inv(M)
D = M\D;
F = M\F;
G = M\G;
H = M\H;
J = M\J;

%Initial condition q in each cell and its interpol. coeffs ro
%(cells indexed by i, each cell has p+1 nodes)
q   = zeros(p+1, Ncells);

%q(:)   = exp(-x(:).^2);
q0 = IC(x,6); % square jump
q = q0;
rho = V\q;

% PARAMETERS 
u = 0.;
Dcoeff = 0.1;

% Lax-Friedrich switch parameters
tauL = (1/u)*(u+abs(u))/2;
tauR = (1/u)*(u-abs(u))/2;

% time step used in time marching
dx = (nodex(2:Ncells)-nodex(1:Ncells-1))/(p+1)^2;
dt = min(min(min(dx))/abs(u), min(min(dx.^2))/Dcoeff);

% final time to integrate advection equation to
T = 8;                  

% number of time steps
Ntsteps = floor(T/dt); 

% TIME STEPPING

for tstep = 1:Ntsteps
  
  sigma = rho;
  
  for rkstage = 4:-1:1 %p:-1:1
    %dsigmadx   = LEGdgderiv(sigma, nodex, tauL, tauR, D, F, G, H, J);
    dsdx       = LEGdgderiv(sigma, nodex,    1,    1, D, F, G, H, J);
    d2sdx2     = LEGdgderiv(dsdx,  nodex,    1,    1, D, F, G, H, J);
    %sigma      = rho + (dt/rkstage)*(-u*dsigmadx+Dcoeff*d2sdx2);
    sigma      = rho + (dt/rkstage)*(Dcoeff*d2sdx2);
  end
  rho = sigma;
  
  if ( ~mod(tstep, 400)  )
    q = V*rho;
    
    plot(x,q0,x,q);    axis([xmin xmax 1 2]);    pause(0.1);

  end
end
