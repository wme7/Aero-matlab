function [x,w,V]=GaussLegendre(n)
% Parameters
x=zeros(n,1);
w=zeros(n,1);
m=(n+1)/2;
xm=0.0;
xl=1.0;

%% Abscisas and Weights
for i=1:m
  z=cos(pi*(i-0.25)/(n+0.5));
  while 1
    p1=1.0;
    p2=0.0; 
    for j=1:n
      p3=p2;
      p2=p1;
      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
    end
    pp=n*(z*p1-p2)/(z*z-1.0);
    z1=z;
    z=z1-p1/pp;
    if (abs(z-z1)<eps), 
        break
    end
  end
  x(i)=xm-xl*z;
  x(n+1-i)=xm+xl*z;
  w(i)=2.0*xl/((1.0-z*z)*pp*pp);
  w(n+1-i)=w(i);
end

%% Legendre Vandermonde Matrix
% compute the value of the legendre polinomial of degree 'k' por column
% values 'x'

% Polynomial Degree 
k = n-1;

% Evaluate sLegendre fucntion
for l = 0:k         % all Polynomial Degrees up to k
    j = l + 1;      % Dummy index
    for i = 1:n     % all Element Points in x
        P(j,i) = legendreP(x(i),l);
    end
end
V = P';             % Transpose our result

return;