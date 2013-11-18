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
format shorte
option.L0 = 1;
option.maxIt = 50;
option.printlevel = 1;
option.elemType = 'WG';
option.plotflag = 0;
option.maxN = 2e4;

err = afemPoisson3(node,elem,pde,bdFlag,option,HB);
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'k-+','eta');
latexerrtable(err.N,[err.H1 err.eta])