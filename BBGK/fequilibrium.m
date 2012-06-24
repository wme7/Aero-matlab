function y = fequilibrium(x) 
u = 0;             % velocity
t = 4.383850;      % temperature
r = 0.2253353;     % fugacity
y = 1 ./ ( exp( (x-u).^2 ./ t) ./ r);
