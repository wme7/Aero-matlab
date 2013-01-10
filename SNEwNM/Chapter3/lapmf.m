% LAMPF Matrix-free negative Laplacian.
% LU=LAPMF(U) is the five-point Laplacian of U.
%
function lu=lapmf(u)
n2=length(u);
n=sqrt(n2);
h2=(n+1)*(n+1);
%
% transform the linear array into a 2d array
%
uu=zeros(n,n);
uu(:)=u;
%
% compute the discrete negative Laplacian on the 2d array
%
luu=4*uu;
if n>1
luu(1:n-1,:)=luu(1:n-1,:)-uu(2:n,:);
luu(2:n,:)=luu(2:n,:)-uu(1:n-1,:);
luu(:,2:n)=luu(:,2:n)-uu(:,1:n-1);
luu(:,1:n-1)=luu(:,1:n-1)-uu(:,2:n);
end
%
% convert to a linear array, divide by h^2
%
lu=luu(:)*h2;
