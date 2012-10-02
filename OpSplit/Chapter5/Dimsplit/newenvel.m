function e=newenvel(f,u,i,j,ncalled)
% Finds the lower convex envelope of the function f(u) between the points
% u(i) and u(j). e is the index of points in the array u where f(u(k))=the
% envelope evaluated at u(k).


%% A quick solver if f is convex
%  %       e=i:j;
%
%% A slower solver for general f
 if nargin<5,
  ncalled=0;
 else
	 ncalled=ncalled+1;
 end;
if i>j,
  warning(14,' Left larger left than right in newenvel');
  k=i; i=j; j=k;
end;
 if i==j,
   e=i; 						%Terminate with a 'shock'
elseif (i+1)==j,
   e=[i j]; 					%Terminate with a 'rarefaction'
 else
   du=u(i)-u(i+1:j);
   df=f(i)-f(i+1:j);
   s=df./du;
   [m k]=min(s);		                %New point on envelope
   e=[i newenvel(f,u,i+k,j,ncalled)];	%Find the rest of the envelope
end;	
