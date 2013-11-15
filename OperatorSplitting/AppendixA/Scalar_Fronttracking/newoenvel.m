function e=newoenvel(f,u,j,i);
%
% Finds the upper envelope of a 'function' f(u) between u(i) and u(j)
%
if i>j,
  [j i]
  error('i>j');
end;
if i==j,
  e=[i]; 	    	%Terminate with a 'shock'
elseif i+1==j,
  e=[i j]; 		%Terminate with a 'rarefaction'
else
  du=u(i)-u(i+1:j);
  df=f(i)-f(i+1:j);
  s=df./du;
  [m k]=max(s);			%New point on envelope
  e=[i newoenvel(f,u,j,i+k)];	%Find the rest of the envelope
end;	
