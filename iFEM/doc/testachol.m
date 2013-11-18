b = ones(size(A,1),1);
tol = 1e-6; maxit = 100;
% tic; L1 = chol(A,'lower'); toc;
% fprintf('\n Incomplete chol decomposition');
% tic; 
% L1 = ichol(A); 
% toc;
% fprintf('\n ichol as Preconditioner');
% tic;
% [x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,@(x)icholpre(x,A,L1,L1'));
% toc;
% fprintf('#dof: %8.0u,  iter: %2.0u\n',size(A,1), it1)
% semilogy(0:it1,rv1./norm(b),'r.');
% 
fprintf('\n Approximate chol decomposition');
tic;
[L2,p,Ac] = achol(A); 
toc;
fprintf('\n Achol as Preconditioner');
tic;
Ap = A(p,p);
[x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,@(r)acholpre(r,Ap,L2,L2',p,Ac));
toc;
fprintf('#dof: %8.0u,  iter: %2.0u\n',size(A,1), it2)
semilogy(0:it2,rv2./norm(b),'b.');
hold on