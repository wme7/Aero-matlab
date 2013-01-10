function w=pde(u)
% PDE  W = PDE(U) is the residual for the convection-diffusion
%      problem in Chapter 3
%
% 
global rhsf;
%
% C=20
%
v=20*u.*(dxmf(u)+dymf(u));
%
% uncomment this line for the unpreconditioned problem
%
%w=-lapmf(u)+20*u.*(dxmf(u)+dymf(u))-rhsf;
%
% uncomment this line for the preconditioned problem
%
w=u+fish2d(v)-rhsf;
