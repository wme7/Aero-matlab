clear all
n = 300; A=randn(n); A=A+A';

fprintf('\nSymmetric eigenvalue problem (size = %d)\n',n);

fprintf('     \n');
fprintf('     n=%d;\n',n);
fprintf('     A=randn(n); A=A+A'';\n');
fprintf('     \n');
fprintf('     [V,D]=syev(A,''mrrr'');\n');
fprintf('     [V,D]=syev(A,''dc'');\n');
fprintf('     [V,D]=syev(A,''bx'');\n');
fprintf('     [V,D]=syev(A,''qr'');\n');
fprintf('     \n');

tic; [V,D]=syev(A,'mrrr'); time=toc;
fprintf('SYEV MRRR  :: time =%4.2f\tnorm(A-V*diag(diag(D))*V'')/norm(A)=%4.2e\tnorm(eye(n)-V''*V)=%4.2e\n',time,norm(A-V*diag(diag(D))*V')/norm(A),norm(eye(n)-V'*V));
tic; [V,D]=syev(A,'dc');   time=toc;
fprintf('SYEV D&C   :: time =%4.2f\tnorm(A-V*diag(diag(D))*V'')/norm(A)=%4.2e\tnorm(eye(n)-V''*V)=%4.2e\n',time,norm(A-V*diag(diag(D))*V')/norm(A),norm(eye(n)-V'*V));
tic; [V,D]=syev(A,'bx');   time=toc;
fprintf('SYEV BX    :: time =%4.2f\tnorm(A-V*diag(diag(D))*V'')/norm(A)=%4.2e\tnorm(eye(n)-V''*V)=%4.2e\n',time,norm(A-V*diag(diag(D))*V')/norm(A),norm(eye(n)-V'*V));
tic; [V,D]=syev(A,'qr');   time=toc;
fprintf('SYEV QR    :: time =%4.2f\tnorm(A-V*diag(diag(D))*V'')/norm(A)=%4.2e\tnorm(eye(n)-V''*V)=%4.2e\n',time,norm(A-V*diag(diag(D))*V')/norm(A),norm(eye(n)-V'*V));
tic; [V,D]=eig(A);   time=toc;
fprintf('MAtlab EIG :: time =%4.2f\tnorm(A-V*diag(diag(D))*V'')/norm(A)=%4.2e\tnorm(eye(n)-V''*V)=%4.2e\n',time,norm(A-V*diag(diag(D))*V')/norm(A),norm(eye(n)-V'*V));
fprintf('     \n');

