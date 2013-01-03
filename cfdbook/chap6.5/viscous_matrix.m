%VISCOUS_MATRIX
% Vectorized matrix generation by evaluation of stencil coefficients
% Output: Bu, Bv
%       |   b5    |
% [B] = |b2 b3  b4|
%       |   b1    |

% Indexing convention in staggered grid:
%  +--k+1--+ 
%  |       |
%  j  jk  j+1
%  |       |
%  +---k---+ 

tic
r = 1/Re;

% Mask for old u, v in Dirichlet boundary points in time-stepping:
masku = ones(size(XU)); maskv = ones(size(XV)); 

%.................Diagonals of viscous matrix Bu for u........................
au1 = zeros(K,J+1); au2=au1; au3=au1; au4=au1; au5=au1;	   % Diagonals of Bu
au1(2:K,:)   = -2*r.*DXU(2:K,:)./(DYU(1:K-1,:) + DYU(2:K,:));
au2(:,2:J+1) =  -r*DY./DX;
au4(:,1:J) = au2(:,2:J+1);  	au5(1:K-1,:) = au1(2:K,:);
au3          = - au1 - au2 - au4 - au5;

%.................Diagonals of viscous matrix Bv for v........................
av1 = zeros(K+1,J); av2=av1; av3=av1; av4=av1; av5=av1;	   % Diagonals of Bv
av1(2:K+1,:)   = -r*DX./DY;
av2(:,2:J)   = -2*r*DYV(:,2:J)./(DXV(:,1:J-1)+DXV(:,2:J));
av4(:,1:J-1) = av2(:,2:J);	av5(1:K,:) = av1(2:K+1,:);
av3          = - av1 - av2 - av4 - av5;

% .................Boundary corrections...............................
jseg = [0,cumsum(nx)];		% j-indices of horizontal segment boundaries
for seg = 1:length(xbc(1,:))
  if (xbc(1,seg) == 1)|(xbc(1,seg) == 2)  % No-slip or inflow at lower boundary
    j = 1+jseg(seg);            au3(1,j) = au3(1,j) + r*DX(1,j)/DY(1,j);
    j = 2+jseg(seg):jseg(seg+1);au3(1,j) = au3(1,j) + 2*r*DXU(1,j)./DYU(1,j);
    j = 1+jseg(seg+1);          au3(1,j) = au3(1,j) + r*DX(1,j-1)/DY(1,j-1);

  elseif xbc(1,seg) == 3	 	  % Outflow at lower boundary
    % Do nothing
  else
    error('Wrong xbc')
  end
  if (xbc(2,seg) == 1)|(xbc(2,seg) == 2)  % No-slip or inflow at upper boundary
    j = 1+jseg(seg);		au3(K,j) = au3(K,j) + r*DX(K,j)/DY(K,j);
    j = 2+jseg(seg):jseg(seg+1);au3(K,j) = au3(K,j) + 2*r*DXU(K,j)./DYU(K,j);
    j = 1+jseg(seg+1);		au3(K,j) = au3(K,j) + r*DX(K,j-1)/DY(K,j-1);          

  elseif xbc(2,seg) == 3		  % Outflow at upper boundary
    % Do nothing
  else
    error('Wrong xbc')
  end
end
kseg = [0,cumsum(ny)];		% k-indices of vertical segment boundaries
for seg = 1:length(ybc(1,:))
  if (ybc(1,seg) == 1)|(ybc(1,seg) == 2)  % No-slip or inflow at left boundary
    k = 1+kseg(seg):kseg(seg+1); 
    au3(k,1) = 0;  masku(k,1) = 0; 
    au1(k,1) = 0; au2(k,1) = 0; au4(k,1) = 0; au5(k,1) = 0;
    k = 1+kseg(seg);            av3(k,1) = av3(k,1) + r*DY(k,1)/DX(k,1);
    k = 2+kseg(seg):kseg(seg+1);av3(k,1) = av3(k,1) + 2*r*DYV(k,1)./DXV(k,1);
    k = 1+kseg(seg+1);          av3(k,1) = av3(k,1) + r*DY(k-1,1)/DX(k-1,1);
             
  elseif ybc(1,seg) == 3	 	  % Outflow at left boundary
    % Do nothing
  else
    error('Wrong ybc')
  end

  if (ybc(2,seg) == 1)|(ybc(2,seg) == 2)  % No-slip or inflow at right boundary
    k = 1+kseg(seg):kseg(seg+1); 
    au3(k,J+1) = 0;  masku(k,J+1) = 0; 
    au1(k,J+1) = 0; au2(k,J+1) = 0; au4(k,J+1) = 0; au5(k,J+1) = 0;
    k = 1+kseg(seg);            av3(k,J) = av3(k,J) + r*DY(k,J)/DX(k,J);
    k = 2+kseg(seg):kseg(seg+1);av3(k,J) = av3(k,J) + 2*r*DYV(k,J)./DXV(k,J);
    k = 1+kseg(seg+1);          av3(k,J) = av3(k,J) + r*DY(k-1,J)/DX(k-1,J);
             
  elseif ybc(2,seg) == 3	 	  % Outflow at left boundary
    % Do nothing
  else
    error('Wrong ybc')
  end
end
for seg = 1:length(xbc(1,:))
  if (xbc(1,seg) == 1)|(xbc(1,seg) == 2)  % No-slip or inflow at lower boundary
    j = 1+jseg(seg):jseg(seg+1); 
    av3(1,j) = 0; maskv(1,j) = 0;
    av1(1,j) = 0; av2(1,j) = 0; av4(1,j) = 0; av5(1,j) = 0;         

  elseif xbc(1,seg) == 3	 	  % Outflow at lower boundary
    % Do nothing
  else
    error('Wrong xbc')
  end
  if (xbc(2,seg) == 1)|(xbc(2,seg) == 2)  % No-slip or inflow at upper boundary
    j = 1+jseg(seg):jseg(seg+1);
    av3(K+1,j) = 0; maskv(K+1,j) = 0; 
    av1(K+1,j) = 0; av2(K+1,j) = 0; av4(K+1,j) = 0; av5(K+1,j) = 0;         

  elseif xbc(2,seg) == 3		  % Outflow at upper boundary
    % Do nothing
  else
    error('Wrong xbc')
  end
end
  
au1 = au1'; au1 = au1(:);	au2 = au2'; au2 = au2(:);
au3 = au3'; au3 = au3(:);	au4 = au4'; au4 = au4(:);
au5 = au5'; au5 = au5(:);
[au1  au2  au3  au4  au5];			% Display stencil on screen
n = (J+1)*K;	
au1 = [au1(J+2:n); zeros(J+1,1)];		% Shifts to accomodate spdiags
au2 = [au2(2:n); 0]; au4 = [0; au4(1:n-1)];  au5 = [zeros(J+1,1); au5(1:n-J-1)];
d = [-J-1;  -1;  0;  1;  J+1];
Bu = invvolu*spdiags([au1  au2  au3  au4  au5], d, n, n);

av1 = av1'; av1 = av1(:);	av2 = av2'; av2 = av2(:);
av3 = av3'; av3 = av3(:);	av4 = av4'; av4 = av4(:);
av5 = av5'; av5 = av5(:);
%[av1  av2  av3  av4  av5]			% Display stencil on screen	
n = J*(K+1);
av1 = [av1(J+1:n); zeros(J,1)];		% Shifts to accomodate spdiags
av2 = [av2(2:n); 0];	av4 = [0; av4(1:n-1)];	av5 = [zeros(J,1); av5(1:n-J)];
d = [-J;  -1;  0;  1;  J];
Bv = invvolv*spdiags([av1  av2  av3  av4  av5], d, n, n);
masku = masku'; masku = masku(:); 
Masku = spdiags(masku,0,length(masku),length(masku));
maskv = maskv'; maskv = maskv(:); 
Maskv = spdiags(maskv,0,length(maskv),length(maskv));

clear au1 au2 au3 au4 au5 av1 av2 av3 av4 av5 d masku maskv
tijd = toc; disp(['viscous_matrix time = ',num2str(tijd)])
