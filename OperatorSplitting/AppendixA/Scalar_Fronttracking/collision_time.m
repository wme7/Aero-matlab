function t=collision_time(s,t,x,T);
% finds the collsion time between x(1) at t(1) with speed s(s) and x(2)
% at t(2) with speed s(2).
%
if nargin<4
  T=1e6;
end;
sl=s(1); sr=s(2);
tl=t(1); tr=t(2);
xl=x(1); xr=x(2);
if sr > sl, 
  t=T+1;
elseif sr < sl, 
  h1=sl*tl;
  h2=sr*tr;
  h3=xl-xr;
  h4=sl-sr;
  h1=(h1-h2-h3);
  t=h1/h4;	
elseif xl==xr&tl==tr,
  t=tl;
else
  t=T+1;
end;
			
