%RIGHT_HAND_SIDE
% Generates contribution of boundary conditions to right-hand side ru, rv for 
% momentum equations

% Functions called: ubd, vbd, pbd

% Indexing convention in staggered grid:
%  +--k+1--+ 
%  |       |
%  j  jk  j+1
%  |       |
%  +---k---+ 

if n < 2, tic, end

tm = t + (omega-1)*dt;
ru = zeros(size(XU)); rv = zeros(size(XV));
for seg = 1:length(xbc(1,:))
  if xbc(1,seg) == 1  			% No-slip at lower boundary
	% Do nothing
  elseif xbc(1,seg) == 2		% Inflow at lower boundary
    j = 1+jseg(seg);  
    ru(1,j) = ru(1,j) +...
      ubd(tm,j,1,'lower')*r*DX(1,j)/(DY(1,j)*volu(1,j));
    for j = 2+jseg(seg):jseg(seg+1)
      ru(1,j) = ru(1,j) +... 
        ubd(tm,j,1,'lower')*2*r*DXU(1,j)./(DYU(1,j).*volu(1,j));
    end
    j = 1+jseg(seg+1);
    ru(1,j) = ru(1,j) +...
      ubd(tm,j,1,'lower')*r*DX(1,j-1)/(DY(1,j-1)*volu(1,j));
  elseif xbc(1,seg) == 3	 	% Outflow at lower boundary
    % Do nothing
  else
    error('Wrong value for xbc in right_hand_side.m, line a')  	
  end
  if xbc(2,seg) == 1  			% No-slip at upper boundary
	% Do nothing
  elseif xbc(2,seg) == 2		% Inflow at upper boundary
    j = 1+jseg(seg);  
    ru(K,j) = ru(K,j) +...
      ubd(tm,j,K+1,'upper')*r*DX(K,j)/(DY(K,j)*volu(K,j));
    for j = 2+jseg(seg):jseg(seg+1);
      ru(K,j) = ru(K,j) +... 
        ubd(tm,j,K+1,'upper')*2*r*DXU(K,j)./(DYU(K,j).*volu(K,j));
    end
    j = 1+jseg(seg+1);
    ru(K,j) = ru(K,j) +...
      ubd(tm,j,K+1,'upper')*r*DX(K,j-1)/(DY(K,j-1)*volu(K,j));
  elseif xbc(2,seg) == 3	 	% Outflow at upper boundary
	% Do nothing
  else
    error('Wrong value for xbc in right_hand_side.m, line b')  	
  end
end

for seg = 1:length(ybc(1,:))
  if ybc(1,seg) == 1			% No-slip at left boundary
    k = 1+kseg(seg):kseg(seg+1);
    ru(k,1) = 0;  
  elseif ybc(1,seg) == 2  		% Inflow at left boundary         
    for k = 1+kseg(seg):kseg(seg+1)
      ru(k,1) = ubd(tm,1,k,'left')/dt; 
    end
    k = 1+kseg(seg);
    rv(k,1) = rv(k,1) +...
      vbd(tm,1,k,'left')*r*DY(k,1)/(DX(k,1)*volv(k,1));
    for k = 2+kseg(seg):kseg(seg+1)
      rv(k,1) = rv(k,1) +...
        vbd(tm,1,k,'left')*2*r*DYV(k,1)./(DXV(k,1).*volv(k,1));
    end
    k = 1+kseg(seg+1);
    rv(k,1) = rv(k,1) +...
      vbd(tm,1,k,'left')*r*DY(k-1,1)/(DX(k-1,1)*volv(k,1));
  elseif ybc(1,seg) == 3	 	% Outflow at left boundary
    for k = 1+kseg(seg):kseg(seg+1)
      ru(k,1) =  ru(k,1) + pbd(tm,xu(1),yu(k),'left')*2./DX(k,1);
    end 
  else
    error('Wrong value for ybc in right_hand_side.m, line c')
  end

  if ybc(2,seg) == 1			% No-slip  at right boundary
    for k = 1+kseg(seg):kseg(seg+1) 
      ru(k,J+1) = 0; 
    end
  elseif ybc(2,seg) == 2  		% Inflow at right boundary
    for k = 1+kseg(seg):kseg(seg+1) 
      ru(k,J+1) = ubd(tm,J+1,k,'right')/dt;
    end
    k = 1+kseg(seg);
    rv(k,J) = rv(k,J) +...
      vbd(tm,J+1,k,'right')*r*DY(k,J)/(DX(k,J)*volv(k,J));
    for k = 2+kseg(seg):kseg(seg+1)
      rv(k,J) = rv(k,J) +...
        vbd(tm,J+1,k,'right')*2*r*DYV(k,J)./(DXV(k,J).*volv(k,J));
    end
    k = 1+kseg(seg+1);          
    rv(k,J) = rv(k,J) +...
      vbd(tm,J+1,k,'right')*r*DY(k-1,J)/(DX(k-1,J)*volv(k,J));             
  elseif ybc(2,seg) == 3	 	  % Outflow at right boundary
    for k = 1+kseg(seg):kseg(seg+1)
      ru(k,J+1) = ru(k,J+1) - pbd(tm,xu(J+1),yu(k),'right')*2./DX(k,J);
    end     
  else
    error('Wrong value for ybc in right_hand_side.m, line d')
  end
end
  
for seg = 1:length(xbc(1,:))
  if xbc(1,seg) == 1  			% No-slip at lower boundary
    j = 1+jseg(seg):jseg(seg+1); 
    rv(1,j) = 0; 
  elseif xbc(1,seg) == 2		% Inflow at lower boundary  
    for j = 1+jseg(seg):jseg(seg+1) 
      rv(1,j) = vbd(tm,j,1,'lower')'/dt;
    end
  elseif xbc(1,seg) == 3	 	% Outflow at lower boundary
    for j = 1+jseg(seg):jseg(seg+1) 
      rv(1,j) = rv(1,j) + pbd(tm,xv(j),yv(1),'lower')*2./DY(1,j);
    end
  end
  if xbc(2,seg) == 1  			% No-slip at upper boundary
    for j = 1+jseg(seg):jseg(seg+1) 
      rv(K+1,j) = 0; 
    end
  elseif xbc(2,seg) == 2		% Inflow at upper boundary
    for j = 1+jseg(seg):jseg(seg+1) 
      rv(K+1,j) = vbd(tm,j,K+1,'upper')/dt;
    end
  elseif xbc(2,seg) == 3	 	% Outflow at upper boundary
    for j = 1+jseg(seg):jseg(seg+1) 
      rv(K+1,j) = rv(K+1,j) - pbd(tm,xv(j),yv(K+1),'upper')*2./DY(K,j);
    end	
  else
    error('Wrong ybc')
  end
end
if n < 2, tijd = toc; disp(['right_hand_side time = ',num2str(tijd)]), end
