function u=testrot(u0,lambda)
  N=length(u0);
  shift=floor(lambda);
  alpha=lambda-shift;
  il=shift+1;
  ir=N-shift-1;
  uh=[u0(1)*ones(1,1:il) u0 u0(N)*ones(1,ir:N)];
  u=alpha*uh(il:il+N-1)+(1-alpha)*uh(il+1:il+N);
  
  