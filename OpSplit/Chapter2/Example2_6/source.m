function h=source(t,u)
scale=1;
h=-scale*u.*(1-u).*(0.5-u)+0.25*u;