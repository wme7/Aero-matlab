%% RATE OF CONVERGENCE OF ADAPTIVE FINITE ELEMENT METHOD USING WG ELEMENT
%
% This example is to show the rate of convergence of the lowest order
% finite element approximation of the second order elliptic equation.

%% Lshape problem
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
% elem = fixorder3(node,elem);
N0 = size(node,1);
HB = zeros(N0,4);
HB(1:N0,1:3) = repmat((1:N0)',1,3); 

bdFlag = setboundary3(node,elem,'Dirichlet');
pde = Lshapedata3;
option.L0 = 1;
option.maxIt = 20;
option.printlevel = 1;
option.elemType = 'WG';
option.plotflag = 1;
option.maxN = 2e4;

[err,time,solver,eqn,node,elem] = afemPoisson3(node,elem,pde,bdFlag,option,HB);
figure;
subplot(1,2,1);
showboundary3(node,elem);
subplot(1,2,2); 
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'k-+','eta');
% latexerrtable(err.N,[err.H1 err.eta])