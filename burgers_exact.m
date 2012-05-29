function burgers_exact
% Plot the exact solution of Burgers' equation
% u_t + u*u_x = 0, u(x,0)=f(x)

u = 0.01:0.01:0.999;

xp = @(u,t) u.*t + sqrt((1 - u)./u) ;
xm = @(u,t) u.*t - sqrt((1 - u)./u) ;

plot(xp(u,0),u,'-b',xm(u,0),u,'-b','LineWidth',2), grid on
xlabel 'x', ylabel 'u'
title ' Solution of u_{t} + uu_{x} = 0'

for t = 0:1:10
plot(xp(u,t),u,'-k',xm(u,t),u,'-k','LineWidth',2), grid on
pause(0.2);
end
