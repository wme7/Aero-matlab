function e=newoenvel(f,u,j,i)
% Finds the lower convex envelope of the function f(u) between the points
% u(i) and u(j). e is the index of points in the array u where f(u(k))=the
% envelope evaluated at u(k).

%% A quick solver if f is convex
%  %       e=[i j];
%
%% A slower solver for general f
if i>j,
 	warning(14,' Larger right than left in newoenvel');
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
 	[m k]=max(s);				%New point on envelope
 	e=[i newoenvel(f,u,j,i+k)];	%Find the rest of the envelope
 end;	
