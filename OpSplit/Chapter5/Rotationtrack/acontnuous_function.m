function z=acontnuous_function(x,y)
r1=(x-0.4).^2+y.^2;
r2=(x+0.3).^2+0.5*y.^2;
z=(sin(pi*min(r1/0.6,1)))+exp(-10*r2);

