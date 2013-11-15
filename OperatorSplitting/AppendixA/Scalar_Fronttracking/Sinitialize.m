function [u,s,t_coll,x]=Sinitialize(u0,x0,T,uflux,fflux);
% Sets up the initial collision times and fronts for front tracking.  
[nc n]=size(u0);
[ur s]=SRiemannsol(u0(1),u0(2),uflux,fflux);
u=[u0(1) ur u0(2)];
[nc nlast]=size(u);
x=x0(1)*ones(size(s));
t_coll=(T+1)*ones([1,nlast]);
for i=2:n-1
  [ur sri]=SRiemannsol(u(nlast),u0(i+1),uflux,fflux);
  u=[u ur u0(i+1)];
  s=[s sri];
  t_coll(nlast)=collision_time([s(nlast-1) sri(1)],[0,0],[x0(i-1),x0(i)]);
  ns=length(sri);
  t_coll=[t_coll (T+1)*ones([1,ns])];
  x=[x x0(i)*ones([1,ns])]; 
  [nc nlast]=size(u);
end;
	
	
