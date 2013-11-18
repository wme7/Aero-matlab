function squareStokesMini
%% 
close all; clear all; clc;
%---------------------- Parameters ----------------------------------------
maxIt = 4; N = zeros(maxIt,1); 
erru = zeros(maxIt,1);  errp = zeros(maxIt,1);
%---------------------- Generate initial mesh -----------------------------
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
for k = 1:3
    [node,elem] = uniformbisect(node,elem);
%     [node,elem] = uniformrefine(node,elem);
end

%% pde
pde = Stokesdata2;
%% 
% option.solver = 'mg';
option.solver = 'direct';
%---------------------- Finite Element Method -----------------------------
for k = 1:maxIt
    [node,elem] = uniformbisect(node,elem);
%     [node,elem] = uniformrefine(node,elem);
    [u,p,A] = StokesMini(node,elem,pde,[],option);
    N(k) = 3*size(node,1);
    uI = pde.exactu(node);
    erru(k) = sqrt((u-uI(:))'*A*(u-uI(:)));
    errp(k) = getL2error(node,elem,@pde.exactp,p);
end
% profile viewer;
%% Plot convergence rates
figure(2);
r1 = showrate(N,erru,2,'-*');
r2 = showrate(N,errp,2,'m-+');
legend('||Du-Du_h||',['N^{' num2str(r1) '}'], ...
       '||p-p_h||',['N^{' num2str(r2) '}'], ...
       'LOCATION','Best');
end 