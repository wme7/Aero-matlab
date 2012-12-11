%
% EFG1D - ONE DIMENSIONAL Galerkin-based MESHLESS PROGRAM FOR SOLVING A 1D BAR
%         SUBJECTED TO A LINEAR BODY FORCE OF MAGNITUDE X WHOSE EXACT SOLUTION IS GIVEN BY
%         u = frac{50}{3} (x - x^3)/E, sigma = frac{50}{3}*(1-3x^2)
%
%         BACKGROUND CELL QUADRATURE IS EMPLOYED TO EVALUATE INTEGRALS
%             - CELLS ARE COINCIDE WITH THE INTERVALS BETWEEN THE NODES
%
%         LAGRANGIAN MULTIPLIER METHOD IS EMPLOYER TO IMPOSE THE ESSENTIAL BOUNDARY CONDITIONS
clear all
% SET UP GLOBAL CONTROL PARAMETERS
scale = 2.4;  % Scale used to determine the radius of support for nodes
nint = 3;     % Order of Gauss quadrature
dx   = 0.1;   % Distance between adjacent nodes
base = 2;     % Basis type - 1: Constant basis;  2: Linear basis; 3: Quadratic basis
type = 2;     % Type of quadrature
%   1: Gauss quadrature;   2: Nodal quadrature;   3: Particle quadrature
WeightType = 'SPLIN';   % Type of weight function ('GAUSS', 'QUART', 'SPLIN','CSRBF')
% SET UP NODAL COORDINATES ALONG BAR, DETERMINE NUMBER OF CELLS
L  = 1.0;              % Length of the bar
xi = [0.0 : dx : L];   % Nodal coordinates
nnodes = length(xi);
ncells = nnodes-1;
% SET MATERIAL PROPERITES
E = 1.0;     % Elastic modulus
area = 1.0;  % Area of cross section
% DETERMINE RADIUS OF SUPPORTS FOR EACH NODE
dm = scale*dx*ones(1,nnodes);
% INITIALIZE MATRICES
K = zeros(nnodes);
P = zeros(nnodes,1);
G = zeros(nnodes,2);
% -----------------------------  Quadrature ----------------------------------------
if type == 1    %  Gauss quadrature
    % LOOP OVER CELLS
    for i = 1:ncells
        
        x1 = xi(i);            % Left point of cell i
        x2 = xi(i+1);          % Right point of cell i
        
        x0 = (x1+x2)/2;  % Coordinate of the mid-point of cell i
        h  = x2-x1;      % Length of cell i
        jac = h/2;       % Jacobian for cell i
        
        [r,w] = Gauss(nint);   % Natural coordinates of Gauss quadrature points and weights
        
        % LOOP OVER GAUSS POINTS
        for j = 1:nint
            
            xq = x0 + h*r(j)/2;   % Coordinates of Gauss quadrature points
            % EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES AT GAUSS POINT xg
            [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xq, dm, WeightType, 0.0);
            % ASSEMBLE DISCRETE EQUATIONS
            K = K + (w(j)*E*area*jac)*(DPHI'*DPHI);
            fbody = area*100*xq;
            P = P + (w(j)*fbody*jac)*PHI';
        end
    end
    
elseif type == 2    %  Nodal quadrature
    
    for j = 1:nnodes
        if j == 1
            w = 0.5*(xi(2)-xi(1));
        elseif j == nnodes
            w = 0.5*(xi(nnodes)-xi(nnodes-1));
        else
            w = 0.5*(xi(j+1)-xi(j-1));
        end
        
        % EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES AT GAUSS POINT xg
        [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xi(j), dm, WeightType, 0.0);
        % ASSEMBLE DISCRETE EQUATIONS
        K = K + (w*E*area)*(DPHI'*DPHI);
        fbody = area*100*xi(j);
        P = P + (w*fbody)*PHI';
    end
elseif type == 3    %  Particle quadrature
    
    %  Loop over nodes
    for j = 1:nnodes
        if j == 1
            w = (xi(2)-xi(1))/6;
        elseif j == nnodes
            w = (xi(nnodes)-xi(nnodes-1))/6;
        else
            w = (xi(j+1)-xi(j-1))/6;
        end
        
        % EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES AT GAUSS POINT xg
        [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xi(j), dm, WeightType, 0.0);
        % ASSEMBLE DISCRETE EQUATIONS
        K = K + (w*E*area)*(DPHI'*DPHI);
        fbody = area*100*xi(j);
        P = P + (w*fbody)*PHI';
    end
    % Generate auxilary points
    xa = zeros(ncells,1);
    xa(1) = (xi(2)-xi(1))/2;
    for j = 2:ncells
        xa(j) = xa(j-1) + (xi(j+1)-xi(j-1))/2;
    end
    % Loop over auxilary points
    for j = 1:ncells
        w = 2*(xi(j+1)-xi(j))/3;
        
        % EVALUATE SHAPE FUNCTIONS AND THEIR DERIVATIVES AT GAUSS POINT xg
        [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xa(j), dm, WeightType, 0.0);
        % ASSEMBLE DISCRETE EQUATIONS
        K = K + (w*E*area)*(DPHI'*DPHI);
        fbody = area*100*xa(j);
        P = P + (w*fbody)*PHI';
    end
    
else
    error('Invalid quadrature type !')
end
% ---------------------------------------------------------------
% ENFORCE BOUNDARY CONDITION USING LAGRANGE MULTIPLIERS
% Prescribed displacement boundary at left end
[PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, 0.0, dm, WeightType, 0.0);
G(1:nnodes,1) = -PHI(1:nnodes)';
% Prescribed displacement boundary at right end
[PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, L, dm, WeightType, 0.0);
G(1:nnodes,2) = -PHI(1:nnodes)';
Q = [0 0];
M = [K G; G' zeros(2)];
% SOLVE FOR NODAL PARAMETERS
d  = M(P',Q)';
% ---------------------------------------------------------------
xg = [0.0 : 0.02 : L];   % Coordinates of points at which the results will be outputed
npts = length(xg);
uh = zeros(npts,1);  % Nodal displacements
sh = zeros(npts,1);  % Nodal stress
for j=1:npts
    [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xg(j), dm, WeightType, 0.0);
    uh(j) = PHI * d(1:nnodes);
    sh(j) = E * DPHI * d(1:nnodes);
end
ue = 50/3*(xg - xg.*xg.*xg)/E;  % Exact solution
se = 50/3*(1 - 3*xg.*xg);
% PLOT RESULTS
figure
subplot(1,2,1);  plot(xg, ue, xg, uh);
subplot(1,2,2);  plot(xg, se, xg, sh);
% Output nodal displacements and stresses
fid1 = fopen('G1DBarDis.dat','w');
fid2 = fopen('G1DBarStr.dat','w');
fprintf(fid1,'%10s%10s%10sn', 'x', 'ue','uh');
fprintf(fid2,'%10s%10s%10sn', 'x', 'se','sh');
for j = 1 : npts
    fprintf(fid1,'%10.6f%10.6f%10.6fn', xg(j), ue(j), uh(j));
    fprintf(fid2,'%10.6f%10.6f%10.6fn', xg(j), se(j), sh(j));
end

fclose(fid1);
fclose(fid2);
% EVALUATE RELATIVE ERROR NORMS BY USING GAUSS QUADRATURE
Luh = 0.0;
%Lue = 0.0;
Lsh = 0.0;
%Lse = 0.0;
% LOOP OVER CELLS
for i = 1:ncells
    
    x1 = xi(i);            % Left point of cell i
    x2 = xi(i+1);          % Right point of cell i
    
    x0 = (x1+x2)/2;  % Coordinate of the mid-point of cell i
    h  = x2-x1;      % Length of cell i
    jac = h/2;       % Jacobian for cell i
    
    [r,w] = Gauss(nint);   % Natural coordinates of Gauss quadrature points and weights
    
    % LOOP OVER GAUSS POINTS
    for j = 1:nint
        
        xq = x0 + h*r(j)/2;   % Coordinates of Gauss quadrature points
        [PHI, DPHI, DDPHI] = MLS1DShape(base, nnodes, xi, 1, xq, dm, WeightType, 0.0);
        uhq = PHI * d(1:nnodes);
        shq = E * DPHI * d(1:nnodes);
        ueq = 50/3*(xq - xq.*xq.*xq)/E;  % Exact solution
        seq = 50/3*(1 - 3*xq.*xq);
        Luh = Luh + w(j)*jac*(uhq-ueq)*(uhq-ueq);
        %       Lue = Lue + w(j)*jac*abs(ueq);
        
        Lsh = Lsh + w(j)*jac*(shq-seq)*(shq-seq);
        %       Lse = Lse + w(j)*jac*abs(seq);
    end
end
Luh = sqrt(Luh)
Lsh = sqrt(Lsh)