%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Solving 1D Wave Equation using CPR
%                 (Z.J. Wang, Feb. 2012)
%
% Equation:             du/dt + du/dx = 0
% Domain:               [-1,1]
% Initial condition:    u(x,0)=-six(pi*x)
% Boundary condition:   periodic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Initilization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; 
close all ;

% k: degree of solution polynomial; k+1 is the numer of solution points
k = 5
% nc: number of elements
nc = 5
% cfl: cfl number
cfl = 0.04
% time: the time to stop the simulation
time = 2

% Other parameters
nf = nc+1;
dx = 2/nc;
dt = cfl*dx
steps = time/dt

% Solution points in the standard element [-1,1]
for i=1:k+1
    x0(i)= -cos((i-1)*pi/k);
end

%solution array
qo=zeros(k+1,nc);
q=zeros(k+1,nc);
qe=zeros(k+1,nc);
xc=zeros(k+1,nc);

% initial condition
for i=1:nc
    for j=1:k+1
        xc(j,i) = -1+(i-1)*dx+(x0(j)+1)/2*dx;
        q(j,i)=-sin(pi*xc(j,i));
    end
end
qe=q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Differential, lifting and reconstruction coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[alfa,coefd,recon,weit] = LiftingAll(k,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Main solver loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale: dksai/dx
scale = 2/dx;

% 3-stage Runge Kutta
for n=1:steps
    qo = q;
    
    % 1st stage
    res = ComputeResidualCPR(q,nc,k,alfa,coefd);
    q = qo-dt*res*scale;

    % 2nd Stage
    res = ComputeResidualCPR(q,nc,k,alfa,coefd); 
    q = 0.75*qo+0.25*(q-dt*res*scale);

    % 3rd stage
    res = ComputeResidualCPR(q,nc,k,alfa,coefd); 
    q = (qo+2*(q-dt*res*scale))/3;
    
    %Plot the solution and the exact solution
    plot(xc,q,xc,qe,'+-');
    
    %pause(0.1)
    if rem(n,10) == 0
        drawnow;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Error computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1_cell = 0;
for i=1:nc
    cell_av_exact=0;
    cell_av_num=0;

    for j=1:k+1
        cell_av_exact = cell_av_exact + weit(j)*qe(j,i);
        cell_av_num   = cell_av_num   + weit(j)*q(j,i);
    end
    L1_cell = L1_cell + abs(cell_av_exact-cell_av_num);
end
L1_cell = L1_cell/nc

L1_node = 0;
for i=1:nc
    for j=1:k+1
        L1_node = L1_node + abs(q(j,i)-qe(j,i));
    end
end
L1_node = L1_node/(nc*(k+1))
