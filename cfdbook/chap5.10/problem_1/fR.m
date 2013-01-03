function ye = fR(t,D,nalpha,nbeta,u,bc,h)	% Right boundary value
alpha = nalpha*pi; beta = nbeta*pi;
if bc == 1
  ye = - h*beta*sin(beta*(1 - u*t))...
    - alpha*exp(-D*alpha*alpha*t)*sin(alpha*(1 - u*t));
else
  ye = 0.0;
end
