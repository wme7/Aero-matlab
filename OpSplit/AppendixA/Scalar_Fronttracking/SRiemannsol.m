function [u,s]=SRiemannsol(ul,ur,uflux,fflux);
%
%   Solves the scalar Riemann problem given by (ul,ur)
%   Output is the intermediate states and the speeds of the front.
%  
if (nargin~=4)
  error('Not correct number of inputs in Riemannsol');
end;
ind=(ul-uflux)==0;
[m i]=max(ind);
if m==0,
  error(' Left state not breakepoint');;
end;
ind=(ur-uflux)==0;
[m j]=max(ind);
if m==0,
  error(' Right state not breakepoint');
end;
if i==j,
  s=[]; u=[];
  return;
elseif i<j
  e=newenvel(fflux,uflux,i,j);
else
  e=newoenvel(fflux,uflux,i,j);
  e=fliplr(e);
end;
du=diff(uflux(e));
df=diff(fflux(e));
s=df./du;
% A costly speed check (for now disabled) 
%	ds=diff(s);
%	if any(ds<0),
%         [s ds]
%	  error('not increasing speeds in Riemannsol');
%	end;
n=length(e);
u=uflux(e(2:n-1));

