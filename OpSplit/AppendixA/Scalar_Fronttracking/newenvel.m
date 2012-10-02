function e=newenvel(f,u,i,j);
%
% Finds the lower envelope of a 'function' f(u) between u(i) and u(j)
%
if i>j,
  error;
end;
if i==j,
  e=[i];              %Terminate with a 'shock'
elseif i+1==j,
  e=[i j]; 	      %Terminate with a 'rarefaction'
else
  du=u(i)-u(i+1:j);
  df=f(i)-f(i+1:j);
  s=df./du;
  [m k]=min(s);		      %New point on envelope
  e=[i newenvel(f,u,i+k,j)];  %Find the rest of the envelope
end;	
