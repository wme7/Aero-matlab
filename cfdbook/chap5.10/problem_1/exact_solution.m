function ye = exact_solution(t,x,D,nalpha,nbeta,u,h)
alpha = nalpha*pi; beta = nbeta*pi;
ye = h*cos(beta*(x - u*t)) + exp(-D*alpha*alpha*t)*cos(alpha*(x - u*t));
