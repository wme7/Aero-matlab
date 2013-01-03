%GRAD_AND_DIV
% Generates pressure gradient matrices Pu, Pv and divergence matrix D
% Vectorized matrix generation by evaluation of stencil coefficients
% Output: Pu, Pv, Du, Dv

J= length(dx); K = length(dy);

%.................Pressure gradient matrix Pu ....................

%        	  |    0    |
% Stencil: [Pu] = |p1 p2  0 |
%        	  |    0    |

% Indexing convention in staggered grid:
%  +--k+1--+ 
%  |       |
%  j  jk  j+1
%  |       |
%  +---k---+ 
tic
p2 = 1./DXU(1,:)'; p1 = zeros(size(p2)); p1(1:J) = -p2(2:J+1); 
Puu = spdiags([p1  p2],[-1;0], J+1, J); Pu = kron(speye(K),Puu);

%......................Boundary corrections....................................
for seg = 1:length(ybc(1,:))
  if (ybc(1,seg) == 1)|(ybc(1,seg) == 2)  % No-slip or inflow at left boundary
    k = 1+kseg(seg):kseg(seg+1); Pu(1+(k-1)*(J+1),:) = 0;
  end
  if (ybc(2,seg) == 1)|(ybc(2,seg) == 2)  % No-slip or inflow at right boundary
    k = 1+kseg(seg):kseg(seg+1);  Pu(k*(J+1),:) = 0; 
  end
end
tijd = toc; disp(['Breakdown of grad_and_div time'])
disp(['  Pu time = ',num2str(tijd)])

%.................Pressure gradient matrix Pv ....................

%        	  |     0   |
% Stencil: [Pv] = |0   p2  0|
%        	  |    p1   |

tic
pp1 = zeros(size(XV)); pp2 = pp1;    		% Diagonals of Pv
pp1(2:K+1,:) = -1./DYV(2:K+1,:); 	pp2 = -pp1;
pp2(1,:) = 1./DYV(1,:); 		pp2(K+1,:) = 0;

%......................Boundary corrections....................................
for seg = 1:length(xbc(1,:))
  if (xbc(1,seg) == 1)|(xbc(1,seg) == 2)  % No-slip or inflow at lower boundary
    j = 1+jseg(seg):jseg(seg+1);     	pp1(1,j) = 0; pp2(1,j) = 0;      
  end
  if (xbc(2,seg) == 1)|(xbc(2,seg) == 2)  % No-slip or inflow at upper boundary
    j = 1+jseg(seg):jseg(seg+1);	pp1(K+1,j) = 0; pp2(K+1,j) = 0; 
  end
end
n = J*(K+1);   p1 = reshape(pp1',n,1);    p2 = reshape(pp2',n,1);
p1 = [p1(J+1:n); zeros(J,1)];		% Shift to accomodate spdiags
Pv = spdiags([p1  p2],[-J;0], n,n-J);
tijd = toc; disp(['  Pv time = ',num2str(tijd)])

%.................u divergence matrix Du ....................

%		  |   0     |
% Stencil: [Du] = |0  p1  p2|  
%		  |   0     |

tic
Duu = spdiags([-ones(J,1) ones(J,1)], [0;1], J,J+1);
Du = kron(spdiags(dy',0,K,K),Duu);
tijd = toc; disp(['  Du time = ',num2str(tijd)])

%.................v divergence matrix Dv ....................

%		  |   p2   |
% Stencil: [Dv] = |0  p1  0|  
%		  |   0    |

tic
Dvv = spdiags([-ones(K,1) ones(K,1)], [0;1], K,K+1);
Dv = kron(Dvv,spdiags(dx',0,J,J));
tijd = toc; disp(['  Dv time = ',num2str(tijd)])

clear  pp1 pp2 p1 p2 Puu Duu Dvv 			% Save storage
