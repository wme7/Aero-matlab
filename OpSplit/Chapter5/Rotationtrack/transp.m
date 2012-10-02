function u=transp(u0,lambda,dir,maxshift)
% Transports u0 a distance lambda and averages over grid (which is assumed
% to be equispaced) Boundary conditions are Neumann. It is easy to modify
% for other boundary conditions.
if nargin<3,
  dir=1;
end;
if (dir>2),
  error(strcat('You can easily extend this to',num2str(dir)','D yourself ...'));
end;
% Some things ensuring that you solve along the columns of u. Must
% transpose if u is a row vector or if dir = 2.
transpose=0;
S=size(u0);
N=S(1);
if ((N==1)||(dir==2)),
  u0=u0';
  transpose=1;
end
S=size(u0);
N=S(1); M=S(2);
Sl=size(lambda);
if ((Sl(1)==Sl(2))&&(Sl(1)==1)),
	lambda=lambda*ones(1,M);
elseif (Sl(2)==1)&&(Sl(1)>1),
	lambda=lambda';
end;
Sl=size(lambda);
if (Sl(2)~=M),
	error('The shift must have the right dimension');
end;
if nargin<4,
	maxshift=ceil(max(abs(lambda)))+1;
end;
shift=floor(lambda);
alpha=lambda-shift;
eshift=ones(maxshift,1);
uh=[eshift*u0(1,1:M); u0; eshift*u0(N,1:M)];
e=ones(N,1);
alpha=e*alpha;   % Making alpha a square matrix the size of u0
uh=colshift(uh,shift);
u=alpha.*uh(maxshift:maxshift+N-1,:)+(1-alpha).*uh(maxshift+1:maxshift+N,:);
if transpose,
	u=u';
end;
  
  