function s = Simpson2(a, b)
%% Simpson Rule for a prestablished f(x)
% a: is the lower limit of integration
% b: is the upper limit of integration
% for this excersice f(x)=x^2

h=(b-a)/2;
x(1)=a;
x(2)=a+(b-a)*1/2;
x(3)=b;

f(1)=myfunction(x(1));
f(2)=myfunction(x(2));
f(3)=myfunction(x(3));

s=h/3*sum([1 4 1].*f);


