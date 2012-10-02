function [u,x]=SFrontTrack(u0,x0,T,u_flux,f_flux);
% Doing front tracking on the initial fronts defined by (u0,x0) in the
% time interval [0 T]. The piecewise linear flux function is in (u_flux,
% f_flux). u0 must take values in the set {u_flux}.
% This version plots the fronts as the Riemann problmes are solved.
  
[u s t_coll,x]=Sinitialize(u0,x0,T,u_flux,f_flux); 
% Initializes the collision times. 
t_start=zeros(size(x));
asm=mean(abs(s));
xlower=min([x(1)-asm*T,x(1)]);
xhigher=max([x(length(x))+asm*T,x(length(x))]);
if xlower==xhigher
  xlower=xhigher-0.5;
  xhigher=xlower+1;
end;
axis([xlower xhigher 0 T]);
set(gca,'Box','on');
%axes('XLim',[xlower xhigher],'YLim',[0 T],'Box','on');
hold on;
[t i]=min(t_coll);
while t<T&~isempty(s),  % The main loop until no more collisions ...
  [ncr nu]=size(u);
  nx=length(x);
  if nx~=nu-1,
    [nx nu];
    error('Wrong length of x vs u in FrontTrack');
  end;
  xc=x(i)+s(i)*(t-t_start(i));
  line([x(i) xc],[t_start(i) t]);       % Plotting the lines to the 	
  line([x(i-1) xc],[t_start(i-1) t]);   % collision about to be solved
  [ur sr]=SRiemannsol(u(:,i-1),u(:,i+1),u_flux,f_flux);  % Solving the RP
  [nc nr]=size(ur); 
  ns=length(sr);
  if ns>0,
    if i>2,
      t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                                 [x(i-2),xc],T+1);
      if t_coll(i-1)<=t,     % The collision is not after present time.
        [-1 t_coll(i-1) t]   % Fatal error.
        [s(i-2) sr(1) t_start(i-2) t]
        error('not increasing t_coll : a');
      end;
    end;
    if i<nu-1,
      t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                                 [xc,x(i+1)],T+1);
      if t_coll(i+1)<=t,  % The collision is not after present time.
        [1 t_coll(i+1) t] % Fatal error.
        [sr(ns) s(i+1)]
        error('not increasing t_coll: b ');
      end;
    end;
  else
    if i>2&i<nu-1,
      t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                                 [t_start(i-2),t_start(i+1)],...
                                 [x(i-2),x(i+1)],T+1);
    end;					  
  end;
  u=[u(:,1:i-1) ur u(:,i+1:nu)];    % Managing the lists ...
  s=[s(1:i-2) sr s(i+1:nu-1)];
  hone=ones([1,ns]);
  x=[x(1:i-2) xc*hone x(i+1:nu-1)];
  t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
  t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
  [t i]=min(t_coll);
end;
n=length(x);
for i=1:n  % Drawing the fronts up to final time T
  line([x(i),x(i)+s(i)*(T-t_start(i))],[t_start(i) T]);
end;
x=x+s.*(T-t_start);
hold off;


