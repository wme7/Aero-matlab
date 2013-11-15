function []=draw(f,a,b)

%
% draws f in the intervall [a,b] 
%

x=linspace(a,b,100);
y=feval(f,x);
plot(x,y,'k-');


print -deps p
