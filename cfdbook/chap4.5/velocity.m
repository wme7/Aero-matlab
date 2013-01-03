function y = velocity(x,b)
y1 = ones(length(x),1);
y1 = y1 - b*sin(pi*x);
y=y1;
