%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Subroutine to compute the residual for the wave equation
%  assuming the solution points include the two interfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = ComputeResidualCPR(q,nc,k,alfa,coefd)

nf = nc+1;
res=zeros(k+1,nc);
ql=zeros(nf);
qr=zeros(nf);

% for each face, the left and right state
for i=2:nc
    ql(i)=q(k+1,i-1);
    qr(i)=q(1,i);
end

% Periodic boundary conditions
i=1;
ql(i)=q(k+1,nc);
qr(i)=q(1,i);

i=nf;
ql(i)=q(k+1,nc);
qr(i)=q(1,1);

for i=1:nc
    for j=1:k+1
        res(j,i)=0.;
        % the derivative term
        for m=1:k+1
          res(j,i) = res(j,i)+ coefd(j,m)*q(m,i);
        end
        % lifting term
        res(j,i)=res(j,i)-alfa(j,1)*(ql(i)-qr(i));
    end    
end