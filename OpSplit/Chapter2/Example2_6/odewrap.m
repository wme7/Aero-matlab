function u=odewrap(source,u0,time)
% Wrapper for matlab's "ode45"
sf=str2func(source);
a=ode45(sf,[0 time],u0);
n=length(a.x);
u=a.y(:,n);
u=reshape(u,size(u0));



