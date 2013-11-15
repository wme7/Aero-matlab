function tr = Trapezoidal(a, b)
%% Simpson Rule for a prestablished f(x)
% a: is the lower limit of integration
% b: is the upper limit of integration
% for this excersice f(x)=x^2

h=(b-a)/2;
x(1)=a;
x(2)=b;

f(1)=x(1)^2;
f(2)=x(2)^2;

tr=h/2*sum([1 1].*f);