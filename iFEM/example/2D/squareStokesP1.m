function squareStokesP1
%STOKE Summary of this function goes here
%   Detailed explanation goes here

close all; clear all; clc;
%figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.7,0.75]);
%---------------------- Parameters ----------------------------------------
maxIt = 4; N = zeros(maxIt,1); 
erru = zeros(maxIt,1);  errp = zeros(maxIt,1);
%---------------------- Generate initial mesh -----------------------------
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
for k = 1:2
    [node,elem] = uniformbisect(node,elem);
end
disp('sss')
%% pde
pde.f = @f; pde.g_D = @g_D; pde.p=@exactp;
%---------------------- Finite Element Method -----------------------------
%
for k = 1:maxIt
    [node,elem] = uniformbisect(node,elem);    
     [u,v,~,A] = StokesP1P1(node,elem,[],@f,@g_D,1);
%     [u,p,A,~,nodef] = StokesP1P0(node,elem,pde,[]);
%     [u,p,A,nodef] = mgStokesP1P0(node,elem,level,pde,option);
    N(k) = length(u) + length(p);
    uI = exactu(nodef);
    erru(k) = sqrt((u-uI(:))'*A{J}*(u-uI(:)));
    center = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
    pI = exactp(center);
    pI = pI-pI(end);
    errp(k)=max(abs(p-pI));
%     errp(k) = norm(p - pI)/sqrt(N);
end
%% Plot convergence rates
figure(2);
showrate2(N,erru,1,'-*','||Du-Du_h||',N,errp,1,'m-+','||p-p_h||');
end  
%---------------------- End of function SQUARE -----------------------------

%---------------------- Sub functions called by CUBE ----------------------
function z = f(p) % load data (right hand side function)
x = p(:,1); y = p(:,2);
z(:,1) = 2^10*( (1-6*x+6*x.^2).*(y-3*y.^2+2*y.^3) ...
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
%--------------------------------------------------------------------------
function z = exactp(p)
x = p(:,1); y = p(:,2);
z = -(2^8)*(2-12*x+12*x.^2).*(y.^2-2*y.^3+y.^4);
end
%--------------------------------------------------------------------------
function z = g_D(p) % Dirichlet boundary condition
z = exactu(p);
end