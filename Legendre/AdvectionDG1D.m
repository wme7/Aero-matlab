%% Discontinuous Galerkin using Legendre polynomials
% for the advection equation q_t + a*q_x = 0
% using u(x,0)=sin(x) and periodic boundary conditions
% advancing in time using TVD-RK4 method. 
% Modify by Manuel Diaz, 2012.09.16

clear all; close all;

%% ORDER OF THE SCHEME
% Order of polynomials used for interpolation:
p = 3;

% Chebychev nodes on the interval [-1,1]
xcheby = sin(0.5*pi*linspace(-1,1,p+1)');

%% DIVIDE THE DOMAIN INTO CELLS
% number of cells to divide interval into
Ncells = 10;

% Physical limits of our Domain
xmin = -pi(); xmax =  pi();

% Set number of nodes along the domain
Nnodes = Ncells+1;
nodex  = linspace(xmin,xmax,Nnodes);

% for each cell of the domain fill with p+1 Chebychev nodes
x = zeros(p+1, Ncells);
for i=1:Ncells
    x(:,i) = ((nodex(i+1)-nodex(i))/2)*xcheby + (nodex(i)+nodex(i+1))/2;
end

%% CREATE THE MATRICES NEEDED FOR THE DG SCHEME
% Vandermonde matrix for the Chebychev nodes
V = legendreVDM(xcheby,p);

% Differentiation matrix
%(acts on the vector of coefficients of q
% and gives coefficients of q_x)
Dhat = legendreDiff(p);

% Mass matrix
M = legendreMass(p);

% Stiffness matrix
D = M*Dhat;

% Modify D (given the in LHS of our equation)
% i.e. D <- D(n,m)+ 2*(-1)^(n+m); 0<=m,n<=p

Dmod = zeros(p+1);
for i = 1:p+1
    for j = 1:p+1
        Dmod(i,j) = D(i,j) + (-1)^((i-1)+(j-1));
    end
end

% Surface term matrix (boundaries of our elements)
S = (-1).^[0:p]'* ones(1,p+1);

%% INITIAL CONDITIONS
% IC q in each cell and its interpol. coefficients ro
% (cells indexed by i, each cell has p+1 nodes)
q   = zeros(p+1, Ncells);
%q(:) = exp(-x(:).^2);
q(:) = sin(x);
rho = V\q;

% PARAMETERS FOR THE ADVECTION EQUATION
a = 1.;

% width of each cell
dx = nodex(2:Nnodes)-nodex(1:Nnodes-1);

% time step for time marching
dt = min(min(dx))/((p+1)^2);

% final time to integrate advection equation to
T_end= 4;

% number of time steps
Ntsteps = floor(T_end/dt);

% coefficients
coeff   = (2*(0:p)+1)'*(a./(dx));

%% TIME STEPPING

for tstep = 1:Ntsteps
    
    sigma = rho; % load IC
    
    for rkstage = 4:-1:1 % RK 4
        % stage 1
        vol_terms = -Dmod*sigma;
        
        % flux coming in from left neighbors
        surfterm =  S*sigma;
        
        idthis = 2:Ncells;      % element's initial index
        idneig = 1:Ncells-1;    % element's last index
        
        vol_terms(:, idthis) = ...
            (dt/rkstage)*coeff(:,idthis).*(vol_terms(:,idthis)+surfterm(:,idneig));
        
        vol_terms(:, 1) = (dt/rkstage)*coeff(:,1).*vol_terms(:,1);
        
        sigma    = rho + vol_terms;
    end
    % Periodic BC
    sigma(:,1)= sigma(:,Ncells);
    
    % Update info
    rho = sigma; 
    
    if ( ~mod(tstep, 100)  ) % compute solution when mod() == 0
        q = V*rho;
        % plot(x(:),abs(q(:)-exp(-(x(:)-dt*tstep).^2)));
        plot(x,q);
        pause(0.1);
    end
end
