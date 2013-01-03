function ye = fL(t,D,nalpha,nbeta,u,h)		% Left boundary value
alpha = nalpha*pi; beta = nbeta*pi;
ye = h*cos(beta*u*t) + exp(-D*alpha*alpha*t)*cos(alpha*u*t);
