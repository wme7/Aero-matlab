function squareStokesCR
%STOKE Summary of this function goes here
%   Detailed explanation goes here

close all; clear all;
profile on
%figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.7,0.75]);
%---------------------- Parameters ----------------------------------------
maxIt = 4; N = zeros(maxIt,1); 
erru = zeros(maxIt,1); errp = zeros(maxIt,1); %errp = zeros(maxIt,1);
%---------------------- Generate initial mesh -----------------------------
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
for i = 1:3
    [node,elem] = uniformbisect(node,elem);
end
%% PDE and options
pde.f = @f;
pde.g_D = @g_D;
option.solver = 'direct';

%% Finite Element Method        
for k = 1:maxIt
   [u,p,edge,eqn] = StokesCR(node,elem,pde,[],option);
    N(k) = length(u)+length(p);
    uI = exactu((node(edge(:,1),:)+node(edge(:,2),:))/2);
%     figure(1); showsolution(node,elem,uI(1:size(node,1)));
    erru(k) = sqrt((u-uI(:))'*eqn.A*(u-uI(:)));
    errp(k) = getL2error(node,elem,@exactp,p);
    [node,elem] = uniformbisect(node,elem);    
%     [node,elem] = uniformrefine(node,elem);
end
%% Plot convergence rates
figure(2); clf;
showrate2(N,erru,2,'-*','||Du_I-Du_h||',N,errp,2,'k-+','|| p - p_h||');
profile viewer
end  
%---------------------- End of function SQUARE -----------------------------

%% Sub functions
function z = f(p) % load data (right hand side function)
x = p(:,1); y = p(:,2);
z(:,1) = 2^10*((1-6*x+6*x.^2).*(y-3*y.^2+2*y.^3) ...
       + (x.^2-2*x.^3+x.^4).*(-3+6*y)...
       - (-3+6*x).*(y.^2-2*y.^3+y.^4) );
z(:,2) = -(2^10)*( (-3+6*x).*(y.^2-2*y.^3 ...
       + y.^4)+(x-3*x.^2+2*x.^3).*(1-6*y+6*y.^2)...
       + (1-6*x+6*x.^2).*(y-3*y.^2+2*y.^3) );
end
%--------------------------------------------------------------------------
function z = exactu(p)
x = p(:,1); y = p(:,2);
z(:,1) = -(2^8)*(x.^2-2*x.^3+x.^4).*(2*y-6*y.^2+4*y.^3);
z(:,2) = 2^8*(2*x-6*x.^2+4*x.^3).*(y.^2-2*y.^3+y.^4);
end
% %--------------------------------------------------------------------------
function z = exactp(p)
x = p(:,1); y = p(:,2);
z = -(2^8)*(2-12*x+12*x.^2).*(y.^2-2*y.^3+y.^4);
end
%--------------------------------------------------------------------------
function z = g_D(p) % Dirichlet boundary condition
z = exactu(p);
end