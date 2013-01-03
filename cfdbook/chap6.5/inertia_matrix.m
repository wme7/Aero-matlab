%INERTIA_MATRIX
% Discretization matrices for inertia term and modification of right-hand sides
% Output: Cu, Cv
%       |   c5    |
% [C] = |c2 c3  c4|
%       |   c1    |

% Indexing convention in staggered grid:
%  +--k+1--+ 
%  |       |
%  j  jk  j+1
%  |       |
%  +---k---+ 

if n< 2, tic, end

au1 = zeros(K,J+1); au2=au1; au3=au1; au4=au1; au5=au1;    % Diagonals of Cu
av1 = zeros(K+1,J); av2=av1; av3=av1; av4=av1; av5=av1;    % Diagonals of Cv

if central == 1		% Central scheme for inertia term
  up = reshape(u0,size(XU')); up = up';		% Two-index ordering of u0 
  vp = reshape(v0,size(XV')); vp = vp';		%  and v0 
  upv = [up;up(K,:)];				% upv: old u in cell vertices	
  upv(2:K,:) = (upv(1:K-1,:).*DYU(1:K-1,:) + upv(2:K,:).*DYU(2:K,:))...
  ./(DYU(1:K-1,:) + DYU(2:K,:));
  vpv = [vp,vp(:,J)];				% vpv: old v in cell vertices  
  vpv(:,2:J) = (vpv(:,1:J-1).*DXV(:,1:J-1) + vpv(:,2:J).*DXV(:,2:J))...
  ./(DXV(:,1:J-1) + DXV(:,2:J));
  au1(2:K,:)   = - vpv(2:K,:)./(2*DYU(2:K,:));
  au2(:,2:J+1) = - up(:,1:J)./(2*DXU(:,2:J+1));
  au4(:,1:J)   =   up(:,2:J+1)./(2*DXU(:,1:J));
  au5(1:K-1,:) =   vpv(2:K,:)./(2*DYU(1:K-1,:));
  au3(:,1)     = - up(:,1)./(2*DXU(:,1));			
  au3(:,J+1)   =   up(:,J+1)./(2*DXU(:,J+1));		
  au3(1:K-1,:) =   au3(1:K-1,:) + au5(1:K-1,:);
  au3(2:K,:)   =   au3(2:K,:)   + au1(2:K,:);
% Further contributions to au3(1,:) and au3(K,:) to be added below depending on 
%  boundary conditions

  av1(2:K+1,:) = - vp(1:K,:)./(2*DYV(2:K+1,:));
  av2(:,2:J)   = - upv(:,2:J)./(2*DXV(:,2:J));
  av4(:,1:J-1) =   upv(:,2:J)./(2*DXV(:,1:J-1));
  av5(1:K,:)   =   vp(2:K+1,:)./(2*DYV(1:K,:));
  av3(1,:)     = - vp(1,:)./(2*DYV(1,:));			
  av3(K+1,:)   =   vp(K+1,:)./(2*DYV(K+1,:));		
  av3(:,2:J)   =   av3(:,2:J)   + av2(:,2:J);
  av3(:,1:J-1) =   av3(:,1:J-1) + av4(:,1:J-1);
% Further contributions to av3(:,1) and av3(:,J) to be added below depending on 
%  boundary conditions.

else			% Upwind scheme for inertia term
  up = (u0 + abs(u0))/2; um = (u0 - abs(u0))/2;
  up = reshape(up,size(XU')); up = up'; um = reshape(um,size(XU')); um = um';
  vp = (v0 + abs(v0))/2; vm = (v0 - abs(v0))/2;
  vp = reshape(vp,size(XV')); vp = vp'; vm = reshape(vm,size(XV')); vm = vm';
  vpv = [vp,vp(:,J)];		% Approximation of v+|v|, v-|v|  
  vmv = [vm,vm(:,J)];		%   in cell vertices
  vpv(:,2:J) = (vpv(:,1:J-1).*DXV(:,1:J-1) + vpv(:,2:J).*DXV(:,2:J))...
  ./(DXV(:,1:J-1) + DXV(:,2:J));
  vmv(:,2:J) = (vmv(:,1:J-1).*DXV(:,1:J-1) + vmv(:,2:J).*DXV(:,2:J))...
  ./(DXV(:,1:J-1) + DXV(:,2:J));
  upv = [up;up(K,:)];		% Approximation of u+|u|, u-|u|  
  umv = [um;um(K,:)];		%   in cell vertices
  upv(2:K,:) = (upv(1:K-1,:).*DYU(1:K-1,:) + upv(2:K,:).*DYU(2:K,:))...
  ./(DYU(1:K-1,:) + DYU(2:K,:));
  umv(2:K,:) = (umv(1:K-1,:).*DYU(1:K-1,:) + umv(2:K,:).*DYU(2:K,:))...
  ./(DYU(1:K-1,:) + DYU(2:K,:));

  au1(2:K,:)   = - vpv(2:K,:)./DYU(2:K,:);
  au2(:,2:J+1) = - up(:,1:J)./DXU(:,2:J+1);
  au4(:,1:J)   =   um(:,2:J+1)./DXU(:,1:J);
  au5(1:K-1,:) =   vmv(2:K,:)./DYU(1:K-1,:);
  au3(:,1)     = - um(:,1)./DXU(:,1);
  au3(:,J+1)   =   up(:,J+1)./DXU(:,J+1);
  au3(:,2:J)   =   (up(:,2:J) - um(:,2:J))./DXU(:,2:J);
  au3(1,:)     = au3(1,:)     + (vpv(2,:) -  vmv(1,:))./DYU(1,:);
  au3(K,:)     = au3(K,:)   + (vpv(K+1,:) - vmv(K,:))./DYU(K,:);
  au3(2:K-1,:) = au3(2:K-1,:) + (vpv(3:K,:) - vmv(2:K-1,:))./DYU(2:K-1,:);

  av1(2:K+1,:) = - vp(1:K,:)./DYV(2:K+1,:);
  av2(:,2:J)   = - upv(:,2:J)./DXV(:,2:J);
  av4(:,1:J-1) =   umv(:,2:J)./DXV(:,1:J-1);
  av5(1:K,:)   =   vm(2:K+1,:)./DYV(1:K,:);
  av3(:,1)     =   (upv(:,2) - umv(:,1))./DXV(:,1);
  av3(:,J)     =   (upv(:,J+1) - umv(:,J))./DXV(:,J);
  av3(:,2:J-1) =   (upv(:,3:J) - umv(:,2:J-1))./DXV(:,2:J-1);
  av3(1,:)     =   av3(1,:) - vm(1,:)./DYV(1,:);
  av3(K+1,:)   =   av3(K+1,:) + vp(K+1,:)./DYV(K+1,:);
  av3(2:K,:)   =   av3(2:K,:) + (vp(2:K,:) - vm(2:K,:))./DYV(2:K,:);
end

% .................Boundary corrections...............................

rvm = zeros(size(XV));	% Contribution to right-hand side
tm = t + (omega-1)*dt;
for seg = 1:length(ybc(1,:))
  if (ybc(1,seg) == 1)|(ybc(1,seg) == 2)  % No-slip or inflow at left boundary
    k = 1+kseg(seg);
    rvm(k,1) = rvm(k,1) + 0.5*up(k,1).*vbd(tm,1,k,'left')*DY(k,1);             
    k = 2+kseg(seg):kseg(seg+1); 
    rvm(k,1) = rvm(k,1) + upv(k,1).*vbd(tm,1,k,'left').*DYV(k,1);             
    k = kseg(seg+1)+1; 
    rvm(k,1) = rvm(k,1) + 0.5*up(k-1,1).*vbd(tm,1,k,'left')*DY(k-1,1);         

  elseif ybc(1,seg) == 3		% Outflow at left boundary
    if central == 1
      k = 1+kseg(seg);
      av3(k,1) = av3(k,1) - 0.5*up(k,1)*DY(k,1)/(DYV(k,1)*DXV(k,1));
      k = 2+kseg(seg):kseg(seg+1); 
      av3(k,1) = av3(k,1) - upv(k,1)./DXV(k,1);
      k = kseg(seg+1)+1; 
      av3(k,1) = av3(k,1) - 0.5*up(k-1,1)*DY(k-1,1)/(DYV(k,1)*DXV(k,1));
    else 
      % Do nothing
    end
  else
    error('Wrong ybc')
  end

  if (ybc(2,seg) == 1)|(ybc(2,seg) == 2)  % No-slip or inflow at right boundary
    k = 1+kseg(seg); 
    rvm(k,J) = rvm(k,J) - 0.5*up(k,J+1).*vbd(tm,J,k,'right')*DY(k,J);          
    k = 2+kseg(seg):kseg(seg+1); 
    rvm(k,J) = rvm(k,J) - upv(k,J+1).*vbd(tm,J,k,'right').*DYV(k,J);           
    k = kseg(seg+1)+1; 
    rvm(k,J) = rvm(k,J) - 0.5*up(k-1,J+1).*vbd(tm,J,k,'right')*DY(k-1,J);     
  
  elseif ybc(2,seg) == 3	 	  % Outflow at right boundary
    if central == 1
      k = 1+kseg(seg);
      av3(k,J) = av3(k,J) + 0.5*up(k,J+1)*DY(k,J)/(DYV(k,J)*DXV(k,J));
      k = 2+kseg(seg):kseg(seg+1); 
      av3(k,J) = av3(k,J) + upv(k,J+1)./DXV(k,J);
      k = kseg(seg+1)+1; 
      av3(k,J) = av3(k,J) + 0.5*up(k-1,J+1)*DY(k-1,J)/(DYV(k,J)*DXV(k,J));
    else 
      % Do nothing
    end
  else
    error('Wrong ybc')
  end
end

rum = zeros(size(XU));			  % Contribution to right-hand side
for seg = 1:length(xbc(1,:))
  if (xbc(1,seg) == 1)|(xbc(1,seg) == 2)  % No-slip or inflow at lower boundary
    j = 1+jseg(seg):jseg(seg+1);  av3(1,j) = 0;
    av1(1,j) = 0; av2(1,j) = 0; av4(1,j) = 0; av5(1,j) = 0;         
    j = 1+jseg(seg); 
    rum(1,j) = rum(1,j) + 0.5*vp(1,j).*ubd(tm,j,1,'lower')*DX(1,j); 
    j = 2+jseg(seg):jseg(seg+1); 
    rum(1,j) = rum(1,j) + vpv(1,j).*ubd(tm,j,1,'lower').*DXU(1,j); 
    j = jseg(seg+1)+1; 
    rum(1,j) = rum(1,j) + 0.5*vp(1,j-1).*ubd(tm,j,1,'lower')*DX(1,j-1); 
  elseif xbc(1,seg) == 3	 	  % Outflow at lower boundary
    if central == 1
      j = 1+jseg(seg);
      au3(1,j) =  au3(1,j) - 0.5*vp(1,j)*DX(1,j)/(DXU(1,j)*DYU(1,j));
      j = 2+jseg(seg):jseg(seg+1);
      au3(1,j) =  au3(1,j) - vpv(1,j)./DYU(1,j); 
      j = jseg(seg+1)+1; 
      au3(1,j) =  au3(1,j) - 0.5*vp(1,j-1)*DX(1,j-1)/(DXU(1,j)*DYU(1,j));    
    else 
      % Do nothing
    end 
  else
    error('Wrong xbc')
  end
  if (xbc(2,seg) == 1)|(xbc(2,seg) == 2)  %No-slip or inflow at upper boundary
    j = 1+jseg(seg):jseg(seg+1);
    av3(K+1,j) = 0; 
    av1(K+1,j) = 0; av2(K+1,j) = 0; av4(K+1,j) = 0; av5(K+1,j) = 0;         
    j = 1+jseg(seg);
    rum(K,j) = rum(K,j) - 0.5*vp(K+1,j).*ubd(tm,j,K,'upper')*DX(K,j);
    j = 2+jseg(seg):jseg(seg+1);
    rum(K,j) = rum(K,j) - vpv(K+1,j).*ubd(tm,j,K,'upper').*DXU(K,j);
    j = jseg(seg+1)+1;
    rum(K,j) = rum(K,j) - 0.5*vp(K+1,j-1).*ubd(tm,j,K,'upper')*DX(K,j-1);
  elseif xbc(2,seg) == 3		  % Outflow at upper boundary
    if central == 1
      j = 1+jseg(seg);
      au3(K,j) =  au3(K,j) + 0.5*vp(K+1,j)*DX(K,j)/(DXU(K,j)*DYU(K,j));
      j = 2+jseg(seg):jseg(seg+1);
      au3(K,j) =  au3(K,j) + vpv(K+1,j)./DYU(K,j); 
      j = jseg(seg+1)+1; 
      au3(K,j) =  au3(K,j) + 0.5*vp(K+1,j-1)*DX(K,j-1)/(DXU(K,j)*DYU(K,j));    
    else 
      % Do nothing
    end 
  else
    error('Wrong xbc')
  end
end

for seg = 1:length(ybc(1,:))
  if (ybc(1,seg) == 1)|(ybc(1,seg) == 2)  %No-slip or inflow at left boundary
    k = 1+kseg(seg):kseg(seg+1);  rum(k,1) = 0;  au3(k,1) = 0; 
    au1(k,1) = 0; au2(k,1) = 0; au4(k,1) = 0; au5(k,1) = 0;
  end
  if (ybc(2,seg) == 1)|(ybc(2,seg) == 2)  %No-slip or inflow at right boundary
    k = 1+kseg(seg):kseg(seg+1);  rum(k,J+1) = 0;   au3(k,J+1) = 0; 
    au1(k,J+1) = 0; au2(k,J+1) = 0; au4(k,J+1) = 0; au5(k,J+1) = 0;
  end
end
for seg = 1:length(xbc(1,:))
  if (xbc(1,seg) == 1)|(xbc(1,seg) == 2)  %No-slip or inflow at lower boundary
    j = 1+jseg(seg):jseg(seg+1);  rvm(1,j) = 0; 
  end
  if (xbc(2,seg) == 1)|(xbc(2,seg) == 2)  %No-slip or inflow at upper boundary
    j = 1+jseg(seg):jseg(seg+1); rvm(K,j) = 0;
  end
end

rum = rum./volu; rvm = rvm./volv;	% Contributions to right-hand side
 
au1 = au1'; au1 = au1(:);	au2 = au2'; au2 = au2(:);
au3 = au3'; au3 = au3(:);	au4 = au4'; au4 = au4(:);
au5 = au5'; au5 = au5(:);
%[au1  au2  au3  au4  au5]			% Display stencil on screen
nn = (J+1)*K;	
au1 = [au1(J+2:nn); zeros(J+1,1)];		% Shifts to accomodate spdiags
au2 = [au2(2:nn); 0];au4 = [0; au4(1:nn-1)];au5 = [zeros(J+1,1); au5(1:nn-J-1)];
d = [-J-1;  -1;  0;  1;  J+1];
Cu = spdiags([au1  au2  au3  au4  au5], d, nn, nn);

av1 = av1'; av1 = av1(:);	av2 = av2'; av2 = av2(:);
av3 = av3'; av3 = av3(:);	av4 = av4'; av4 = av4(:);
av5 = av5'; av5 = av5(:);
%[av1  av2  av3  av4  av5]			% Display stencil on screen	
nn = J*(K+1);
av1 = [av1(J+1:nn); zeros(J,1)];		% Shifts to accomodate spdiags
av2 = [av2(2:nn); 0];	av4 = [0; av4(1:nn-1)];	av5 = [zeros(J,1); av5(1:nn-J)];
d = [-J;  -1;  0;  1;  J];
Cv = spdiags([av1  av2  av3  av4  av5], d, nn, nn);

%clear au1 au2 au3 au4 au5 av1 av2 av3 av4 av5 d nn

if n<2, tijd = toc; disp(['inertia_matrix time = ',num2str(tijd)]), end
