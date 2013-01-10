% ATANDEMO  This program creates the ArcTangent example in Chapter 2.
%
% C. T. Kelley, November 13, 2002.
%
x=10;
tol=[1.d-6,1.d-6];
%
% Solve the problem with Newton's method.
%
params=[40, 1, 0,0];
[an, errsn, ierrn]=nsold(x, 'fatan', tol, params);
%
% Set the default Jacobian updating strategy.
%
params=[40, 1000, .5, 0];
[ac, errsd, ierrd]=nsold(x, 'fatan', tol, params);
inewt=length(errsn(:,1)); cnewt=0:inewt-1;
idefault=length(errsd(:,1)); cdefault=0:idefault-1;
%
% Plot the residual histories.
%  
semilogy(cnewt,errsn(:,1),'-',cdefault,errsd(:,1),'--')
legend('Newton','default')
