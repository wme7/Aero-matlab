function [h,hjac]=heq(x)
% HEQ  Chandrasekhar H-equation residual
%      Jacobian uses precomputed data for fast evaluation
%
%      [H, HJAC] = HEQ(X) returns the nonlinear residual
%      H and (optionally) the Jacobian.
%
% Be sure and store the correct data in the global array A_heq.
%
global A_heq;
n=length(x);
h=ones(n,1)-(A_heq*x);
ph=ones(n,1)./h;
h=x-ph;
if nargout==2
    hjac=(ph.*ph)*ones(1,n);
    hjac=A_heq.*hjac;
    hjac=eye(n)-hjac;
end
