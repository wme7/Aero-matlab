clear all
m=40; n = 400; nrhs=1; A=randn(m,n); b=randn(m,nrhs);

fprintf('\nUnderdetermined system of size %d-by-%d, %d right-hand side(s) \n',m,n,nrhs);

fprintf('     \n');
fprintf('     m=40; n = 400; nrhs=1;\n');
fprintf('     A=randn(m,n); b=randn(m,nrhs);\n');
fprintf('     \n');
fprintf('     [x]=gels(A,b) \n');
fprintf('     [x]=gels(A,b,''svd'') \n');
fprintf('     [x]=gels(A,b,''dc'') \n');
fprintf('     [x]=gels(A,b,''y'') \n');
fprintf('     \n');

tic; [x]=gels(A,b); time=toc;
fprintf('     GELS         :: time =%4.2fs\tnorm(x-A''*(A''\\x))=%4.2e\tnorm(A*x-b)=%4.2e\n',time,norm(x-A'*(A'\x)),norm(A*x-b));

tic; [x]=gels(A,b,'svd'); time=toc;
fprintf('     GELS SVD     :: time =%4.2fs\tnorm(x-A''*(A''\\x))=%4.2e\tnorm(A*x-b)=%4.2e\n',time,norm(x-A'*(A'\x)),norm(A*x-b));

tic; [x]=gels(A,b,'dc'); time=toc;
fprintf('     GELS DC      :: time =%4.2fs\tnorm(x-A''*(A''\\x))=%4.2e\tnorm(A*x-b)=%4.2e\n',time,norm(x-A'*(A'\x)),norm(A*x-b));
 
tic; [x]=gels(A,b,'y'); time=toc;
fprintf('     GELS Y       :: time =%4.2fs\tnorm(x-A''*(A''\\x))=%4.2e\tnorm(A*x-b)=%4.2e\n',time,norm(x-A'*(A'\x)),norm(A*x-b));

tic; x=A\b; time=toc;
fprintf('     Matlab \\     :: time =%4.2fs\t                          \tnorm(A*x-b)=%4.2e\n',time,norm(A*x-b));


clear all
m=400; n = 40; nrhs=1; A=randn(m,n); b=randn(m,nrhs);

fprintf('\nOverdetermined system of size %d-by-%d, %d right-hand side(s) \n',m,n,nrhs);

fprintf('     \n');
fprintf('     m=%d; n=%d; nrhs=%d;\n',m,n,nrhs);
fprintf('     A=randn(m,n); b=randn(m,nrhs);\n');
fprintf('     \n');
fprintf('     [x]=gels(A,b) \n');
fprintf('     [x]=gels(A,b,''svd'') \n');
fprintf('     [x]=gels(A,b,''dc'') \n');
fprintf('     [x]=gels(A,b,''y'') \n');
fprintf('     \n');

tic; [x]=gels(A,b); time=toc;
fprintf('     GELS         :: time =%f\tnorm(x-A''*(A''\\x))=%e\tnorm(A''*(A*x-b))=%e\n',time,norm(x-A'*(A'\x)),norm(A'*(A*x-b)));

tic; [x]=gels(A,b,'svd'); time=toc;
fprintf('     GELS SVD     :: time =%f\tnorm(x-A''*(A''\\x))=%e\tnorm(A''*(A*x-b))=%e\n',time,norm(x-A'*(A'\x)),norm(A'*(A*x-b)));

tic; [x]=gels(A,b,'dc'); time=toc;
fprintf('     GELS DC      :: time =%f\tnorm(x-A''*(A''\\x))=%e\tnorm(A''*(A*x-b))=%e\n',time,norm(x-A'*(A'\x)),norm(A'*(A*x-b)));

tic; [x]=gels(A,b,'y'); time=toc;
fprintf('     GELS Y       :: time =%f\tnorm(x-A''*(A''\\x))=%e\tnorm(A''*(A*x-b))=%e\n',time,norm(x-A'*(A'\x)),norm(A'*(A*x-b)));

tic; x=A\b; time=toc;
fprintf('     Matlab \\     :: time =%f\tnorm(x-A''*(A''\\x))=%e\tnorm(A''*(A*x-b))=%e\n',time,norm(x-A'*(A'\x)),norm(A'*(A*x-b)));
fprintf('     \n');
