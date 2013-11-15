function [u,s]=ScalarRiemannsol(ul,ur,delta,f)
%
% Finds the solution of the Riemann problem given by (ul,ur) with a
% parameter delta, and flux 'f'. The results given in 
% (u,s), where length(u)=length(s)-1. This is done by computing the
% lower/upper envelope.
%
if (ul==ur)
    s=[];
    u=[];
    return;
end;
if (abs(ul-ur)<=delta)
    u=[];
    fl=feval(f,ul);
    fr=feval(f,ur);
    s=(fr-fl)/(ur-ul);
    return;
end;
ulow=min(ul,ur);
uhigh=max(ul,ur);
n=ceil((uhigh-ulow)/delta);
uflux=linspace(ulow,uhigh,n); 
fflux=feval(f,uflux);
if ul<ur
	i=1;
	j=n;
else
	i=n;
	j=1;
end;
if i<j
  e=newenvel(fflux,uflux,i,j);
else
  e=newoenvel(fflux,uflux,i,j);
  e=fliplr(e);
end;
du=diff(uflux(e));
df=diff(fflux(e));
s=df./du;
u=uflux(e(2:length(e)-1));

