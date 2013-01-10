function pu=pderight(u)
% PDERIGHT 
% PU = PDERIGHT(U) is the nonlinear residual of the 
% pde example with right preconditioning and C=20.
%
%
global rhsf prhsf;
%
% Apply the preconditioner first.
%
w=fish2d(u);
%
% Low-order nonlinear term. 
%
v=20*w.*(dxmf(w)+dymf(w));
pu=u+v-rhsf;
