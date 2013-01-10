function w=pdeleft(u)
% PDELEFT W=PDELEFT(U) is the nonlinear residual of the left 
%         preconditioned pde example with C=20.
% 
global rhsf prhsf
%
% Compute the low-order nonlinear term.
%
v=20*u.*(dxmf(u)+dymf(u));
%
% Apply fish2d to the entire pde.
%
w=u+fish2d(v)-prhsf;
