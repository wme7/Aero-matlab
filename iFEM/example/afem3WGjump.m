%% JUMPMG2 Diffusion equation with jump coefficients in three dimensions 
%
% The purpose of this subroutine is to test the performance of V-cycle
% multigrid for solving the linear system of algebraic equations arising
% from linear finite element discretization of the elliptic partial
% differential equation with jump coefficients.
%
% $-\nabla \cdot (\omega\nabla u) = f $  in  $\Omega=(-1,1)^3$ 
% $u = 1$ on $x==1$ and $u=0$ on $x==-1$, 
% $\omega\nabla u \cdot n = 0 on other faces.
%
% The diffusion coefficent $\omega$ is piecewise constant with large jump.
%
% jumpmgdata1
%* $\omega(x) = 1/\epsilon$ if $x\in (0, 1)^3$ and 
%* $\omega = 1$ otherwise.  
%
% jumpmgdata2
%* $\omega(x) = 1$ if $x\in (-0.5, 0)^3$ or $x\in (0,0.5)^3$ and 
%* $\omega = \epsilon$ otherwise.  
%
% See also jumpMG1  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; clear all
global epsilon
epsilon = 1e-4;

%% Jump coefficients
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');
pde = jumpmgdata1;

option.L0 = 3;
option.maxIt = 50;
option.printlevel = 1;
% option.elemType = 'WG';
option.plotflag = 1;
option.maxN = 1e4;

[err,time,solver,eqn,node,elem] = afemPoisson3(node,elem,pde,bdFlag,option,HB);

figure;
subplot(1,2,1);
showboundary3(node,elem,'~(x<=0 & y<=0)');
subplot(1,2,2); 
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'k-+','eta');