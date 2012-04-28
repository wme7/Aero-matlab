%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: degree of the solution polynomial
% x0: the solution points in the standard element [-1, 1]
%
% Outputs:
% alfa: lifting coefficients
% coefd: for the derivative du/dksai
% recon: reconstruction coefficients for the left and right interfaces
% weit: quadrature weight (to computate the average)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alfa, coefd, recon, weit] = LiftingAll(k,x0)
alfa=zeros(k+1,2);
coefd=zeros(k+1,k+1);
recon=zeros(2,k+1);
weit=zeros(k+1,1);
M = zeros(k+1,k+1);

% symbolic manipulation
x = sym('x');

% Lagrange polynomials
for i=1:k+1
    l(i)=x/x;
    for j=1:k+1
        if(i ~= j) 
            l(i)=l(i)*(x-x0(j))/(x0(i)-x0(j));
        end
    end
end

% the mass matrix
for i=1:k+1
    weit(i)=int(l(i),-1,1)/2;
    for j=1:i
        M(i,j)=int(l(i)*l(j),-1,1);
        M(j,i)=M(i,j);
    end
end
M=double(M);
weit = double(weit)
%M=vpa(M);
M

for i=1:k+1
  D(i)=diff(l(i),x);
  for j=1:k+1
    coefd(j,i) = subs(D(i),x0(j));
  end
end

for j=1:k+1
  N(j,1)=subs(l(j),-1);  
  N(j,2)=subs(l(j),1);
  recon(1,j)=subs(l(j),-1);
  recon(2,j)=subs(l(j),1);
end
%N
    
alfa=double(inv(M)*N);
