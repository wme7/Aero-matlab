function [fullu,w,u,tn] = CahnHilliardP1(node,elem,pde,w0,u0,tao,times)
tic;
N = size(node,1); NT = size(elem,1);
dquadorder = 5; %the default 
%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);
[lambda,weight] = quadpts(dquadorder);
      phi = lambda;                 
	nQuad = size(lambda,1);
%% Assemble stiffness and fmass matrix and jacobi matrix M0;
A  = sparse(N,N);  
M  = sparse(N,N);
M0 = sparse(N,N);
for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j dxdy$ 
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        
        Mij = 1/12.*area; %according the rule, you can check some finite element books.
                          %when i=j, the answer is differenet . 
          M0ij=0;         
          for p = 1:nQuad 
           M0ij = M0ij + weight(p)*phi(p,i)*phi(p,j)*pde.df2du2((phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
          end
           M0ij = M0ij.*area;                
            if (j==i) 
            A  = A  + sparse(elem(:,i),elem(:,j),Aij,N,N);
            Mij=2*Mij;
            M  = M   + sparse(elem(:,i),elem(:,j),Mij,N,N);
            M0 = M0  + sparse(elem(:,i),elem(:,j),M0ij,N,N);
            else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],N,N);     
            M = M + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Mij; Mij],N,N);
            M0 = M0 + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [M0ij; M0ij],N,N);
            end        
    end
end




  bt = zeros(NT,3);
   for p = 1:nQuad
		% quadrature points in the x-y coordinate
        for i = 1:3
        bt(:,i) = bt(:,i)+weight(p)*phi(p,i)*pde.dfdu((phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
        end
   end
    bt = bt.*repmat(area,1,3);
   rf2 = accumarray(elem(:),bt(:),[N 1]);  % the nonlinear item;
   

clear Aij Mij M0ij Dphi
global efsilonsquare

 
L = [tao*A M; M -efsilonsquare*A];   %the unknow U = [w;u];
 
 rf1   = M*u0; 
 rf    = [rf1;rf2];          % the right hand,including two items;


    U0   = [w0;u0];  % keep the two variables.
   Itmax = 6;        % the max of newton iteration;
     res = L*U0-rf;  % the residual ;
     
 for k  = 1:times
     
    
     tn = tao*k;
%*************************newton method begin******************************
    ii = 0; 
     while(ii<Itmax)
     ii = ii+1;
     Jacobi = L - [sparse(N,N) sparse(N,N);sparse(N,N) M0];
     U0 = U0 - Jacobi\res;
     u0 = U0(N+1:2*N);
     M0 = sparse(N,N);  %Jacobi matrix;
for i = 1:3
    for j = i:3 
         M0ij = zeros(NT,1);
        for p = 1:nQuad
%		quadrature points in the x-y coordinate
        M0ij = M0ij + weight(p)*phi(p,i)*phi(p,j).*pde.df2du2(phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3)));
        end
        M0ij = M0ij.*area;
        if(j==i)
        M0 = M0 + sparse(elem(:,i),elem(:,j),M0ij,N,N);
        else
        M0 = M0 + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                      [M0ij; M0ij],N,N); 
        end
     end
end    


  bt = zeros(NT,3);
   for p = 1:nQuad
		% quadrature points in the x-y coordinate
        for i = 1:3
        bt(:,i) = bt(:,i)+weight(p)*phi(p,i)*pde.dfdu((phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
        end
   end
    bt = bt.*repmat(area,1,3);
   rf2 = accumarray(elem(:),bt(:),[N 1]);
   res = L*U0 -[rf(1:N);rf2];
end
     
 
rf(1:N) = M*u0;



% 
% %************************* newton method end******************************

%************************* direct method begin******************************

% % 
%        bt = zeros(NT,3);
%    for p = 1:nQuad
% 		quadrature points in the x-y coordinate
%      pxy = lambda(p,1)*node(elem(:,1),:) ...
% 			+ lambda(p,2)*node(elem(:,2),:) ...
% 			+ lambda(p,3)*node(elem(:,3),:);
% 		fp = pde.f(pxy,tn);
%         for i = 1:3
%         bt(:,i) = bt(:,i)+weight(p)*phi(p,i)*fp;
%         end        
%    end
%     bt = bt.*repmat(area,1,3);
%    rf2 = accumarray(elem(:),bt(:),[N 1]);
%    
%     rf(1:N) = M*u0+tao*rf2;
%        bt = zeros(NT,3);
%    for p = 1:nQuad
% 		quadrature points in the x-y coordinate
%         s = pde.dfdu((phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
%         s = (phi(p,1)*U0(elem(:,1))+phi(p,2)*U0(elem(:,2))+phi(p,3)*U0(elem(:,3))).^3-...
%         (phi(p,1)*U0(elem(:,1))+phi(p,2)*U0(elem(:,2))+phi(p,3)*U0(elem(:,3)));
%        
%         for i = 1:3
%       bt(:,i) = bt(:,i)+weight(p)*phi(p,i).*s;
%         end
%    end
% 
%            bt = bt.*repmat(area,1,3);
%           rf2 = accumarray(elem(:),bt(:),[N 1]);
%   rf(N+1:2*N) = rf2;
%           
% compute the right term end. 
% 
% U0 = L\rf;
% 
% u0 = U0(N+1:2*N);
% u0
%************************* direct method end*******************************
fullu(:,k+1) = u0;
 k
% tn
% format long e
 end
u = U0(N+1:2*N);
w = U0(1:N);

%  x   = node(:,1);
%  y   = node(:,2);
% %plot3(x,y,u);
% % surf(x,y,u);
% % view(0,90)
% z   = full(u);
% x0  = reshape(x,sqrt(N),sqrt(N));
% y0  = reshape(y,sqrt(N),sqrt(N));
% z0  = reshape(z,sqrt(N),sqrt(N));
%  surf(x0,y0,z0);
% [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),200)',linspace(min(y),max(y),200),'v4');%插值
% figure,contourf(X,Y,Z); %等高线图
% colorbar


end % end of CahnHilliard
