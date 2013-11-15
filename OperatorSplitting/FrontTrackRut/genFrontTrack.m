%% genFrontTrack
function [u,x]=genFrontTrack(u0,x0,T,delta,varargin)
%
%	Does front tracking solution of the initial data 'u0(1:n)', 
%   with discontinuities at the points 'x0(1:n-1) , in the time interval 
%   [T(1) T(2)]. If length(T)=1, then we track in the interval [0,T].
%   using 'delta' as a parameter in the Riemann solver.
%
%     [u,x]=genFrontTrack(u0,x0,T,delta,varargin);
%
%   The following parameters can be set: 
%           'display'  :'lines'      : for drawing the fronts as lines
%	 	    'polygons'   : for drawing colored polygons, must be followed
%                          by the name of the colorfunction.
%
%                          Can be followed by 'fixed_interval', which
%                          again must be followed by a fixed interval
%                          for the colfunc.
%           'none'      : no graphics, default.
%
%           'wbar'      : no graphics, but a waitbar will be displayed
%
%           'save'    : creates output that can be saved, default not.
%
% 			'colors'  : a colormap to be used for colors, must be
%					    followed by the colormap.
%
%           'periodic': if periodic boundary, must be followed by
%                       the interval '[a,b]', where a<x0<b. Also
%                       u0(:,1)=u0(:,length(u0)) in this case.
%
%           'Riemann_solver': specifies the name of the Riemann solver,
%                             followed by the name of the solver
%                             used. This must have the syntax
%                             [u s]=Riemannsol(ul,ur,delta); The default is
%                             'Riemannsol'.
%            'Axis':    followed by a handle to an axis which must already
%                       exist. Must be
%                       set for all options except 'display' = 'none'.
%


%% Preliminaries ....
COLL_TOL=1e-12;  % The tolerance for collision times. 

l=length(varargin);
plain=1; lines=0; polygons=0; savegraphics=0; cmap_set=0; Periodic=0;
Rsolver='Riemannsol'; wb=0; fixedcolinterval=0;
plotAxis=[];
if l==0,
   plain=1;
else
   i=1;
   while i<=l
      if strcmp(varargin(i),'display')
         i=i+1;
         if strcmp(varargin(i),'polygons')
            polygons=1;
            plain=0;
            i=i+1;
            colfunc=varargin{i};
            i=i+1;
            if strcmp(varargin(i),'fixed_interval')
               i=i+1;
               cint=varargin(i);
               pmin=cint{1}(1);
               pmax=cint{1}(2);
               fixedcolinterval=1;
            end;
         elseif strcmp(varargin(i),'lines')
            lines=1;
            plain=0;
         elseif strcmp(varargin(i),'none')
            plain=1;
         elseif strcmp(varargin(i),'wbar')
            wb=1;
            plain=1;
         else
            error('No such option for display');
         end;
      elseif strcmp(varargin(i),'save')
         savegraphics=1;
      elseif strcmp(varargin(i),'colors')
         i=i+1;
         cmap=varargin(i);
         cmap=cmap{1,1};
         cmap_set=1;
      elseif strcmp(varargin(i),'periodic')
         Periodic=1;
         i=i+1;
         bndr=varargin(i);
         xleft=bndr{1}(1);
         xright=bndr{1}(2);
      elseif strcmp(varargin(i),'Riemann_solver')
         i=i+1;
         Rsolver=varargin{i};
      elseif strcmp(varargin(i),'Axis')
         i=i+1;
         plotAxis=varargin{i};
      end;
      i=i+1;
   end;
end;
if isempty(plotAxis)&&(plain==0)
   figure;
   dx0=x0(length(x0))-x0(1);
   axis([x0(1)-dx0  x0(length(x0))+dx0 0 T]);
   plotAxis=gca;
   set(plotAxis,'Box','on','Drawmode','fast');
   display('You should specify an axis in which to plot');
end;
if length(T)==1,
   T=T(1); startTime=0.0;
else
   startTime=T(1);
   T=T(2);
end;

%% Pure front tracking
% To avoid lots of 'if' tests I chose to create rather bulky code
if plain, % READ THIS CASE IF CONCERNED WITH FRONT TRACKING, the other cases differ only in
   % graphics option, such as plotting in various modes.
   
   
   % Initializes the fronts and collision times
   if ~Periodic,
      [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
         'start_time',startTime);
   else
      [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
         'periodic',xleft,xright,...
         'start_time',startTime);
   end;
   % 'u' are the states AFTER the initial Riemann problems have been solved
   % 's' are the speeds of the discontiuity between 'u(i-1)' and 'u(i)'
   % 't_coll' is the collision time of fronts at 'x(i-1)' and 'x(i)'
   t_start=startTime*ones(size(x));
   % The starting times of the fronts
   
   [t i]=min(t_coll);
   % The first collision time is between fronts i-1 and i
   % first=1;
   % Continues tracking while there are more collsions before T, and more
   % discontinuities
   ncoll=0;
   if wb,
      hbar=waitbar(0,'Tracking fronts...');
   end;
   while any(t<T)&&~isempty(s), % Tracking
      ncoll=ncoll+1;
      nu=size(u,2);
      nx=length(x);
      if nx~=nu-1,
         display([nx nu]);
         error('Wrong length of x vs u in FrontTrack');
      end;
      % Finds the collsion point
      xc=x(i)+s(i)*(t-t_start(i));
      
      % solving the Riemann problem, giving intermediate states 'ur'
      % and speeds 'sr'
      if ((Periodic)&&(i==2)),
         tcoll_left=collision_time([0.0 s(3)], [t t_start(3)], [xleft x(3)],T+1);
         tcoll_right=collision_time([s(nu-2) s(2)], [t_start(nu-2) t],...
            [x(nu-2) xright],T+1);
         x=[xleft x(3:nu-2) xright xright];
         s=[0 s(3:nu-2) s(2) 0];
         t_start=[t t_start(3:nu-2) t t];
         u=[u(:,3) u(:,3:nu-1) u(:,3) u(:,3)];
         t_coll=[T tcoll_left t_coll(4:nu-2) tcoll_right T T];
      elseif ((Periodic)&&(i==nu-1)),
         tcoll_right=collision_time([s(nx-2) 0.0],[t_start(nx-2) t],...
            [x(nx-2) xright],T+1);
         tcoll_left=collision_time([s(nx-1),s(2)],[t t_start(2)],...
            [xleft, x(2)],T+1);
         x=[xleft xleft x(2:nx-2) xright];
         t_start=[t t t_start(2:nx-2) t];
         s=[0 s(nx-1) s(2:nx-2) 0];
         u=[u(:,nu-2) u(:,nu-2) u(:,2:nu-2) u(:,nu-2)];
         t_coll=[T,T, tcoll_left,t_coll(3:nu-3),tcoll_right,T];
      else
         [ur sr]=feval(Rsolver,u(:,i-1),u(:,i+1),delta);
         nr = size(ur,2);
         ns = length(sr);
         % updates the collision vector
         if ns>0,
            if i>2,
               t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                  [x(i-2),xc],T+1);
               if t_coll(i-1)<t-COLL_TOL,
                  display([xc t_coll(i-1) t]);
                  display([s(i-2) sr(1) t_start(i-2) t]);
                  error('not increasing t_coll');
               end;
            end;
            if i<nu-1,
               t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                  [xc,x(i+1)],T+1);
               if t_coll(i+1)<t-COLL_TOL,
                  display([xc t_coll(i+1) t]);
                  display([sr(ns) s(i+1)]);
                  error('not increasing t_coll');
               end;
            end;
         else
            ur=[]; sr=[]; nr=0;
            if (i>2)&&i<(nu-1),
               t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                  [t_start(i-2),t_start(i+1)],...
                  [x(i-2),x(i+1)],T+1);
            end;
         end;
         % Updates the speeds
         s=[s(1:i-2) sr s(i+1:nu-1)];
         hone=ones([1,ns]);
         % Updates the locations
         x=[x(1:i-2) xc*hone x(i+1:nu-1)];
         % Updates the starting points
         t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
         % Updates the values and the collsion times
         if ns,
            u=[u(:,1:i-1) ur u(:,i+1:nu)];
            t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
         else
            u=[u(:,1:i-1) u(:,i+2:nu)];
            t_coll=[t_coll(1:i-1) t_coll(i+2:nu)];
         end;
      end;
      % Finds the next collision
      if wb,
         waitbar(t/T,hbar);
      end;
      [t i]=min(t_coll);
   end;
   % Moves the discontinuities up to T
   x=x+s.*(T-t_start);
   if wb,
      close(hbar);
   end;
else
   %% Drawing with simple lines
   if lines,
      if savegraphics,
         %% Saving the graphics, takes some time and memory...
         if ~Periodic,
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime);
         else
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime,'periodic',xleft,xright);
         end;
         t_start=startTime*ones(size(x));
         
         %lims=get(plotAxis,'XLim');
         % xlower=lims(1); xhigher=lims(2);
         hold on;
         [t i]=min(t_coll);
         while (t<T)&&~isempty(s),
            nu=size(u,2);
            nx=length(x);
            if nx~=nu-1,
               display([nx nu]);
               error('Wrong length of x vs u in FrontTrack');
            end;
            xc=x(i)+s(i)*(t-t_start(i));
            line([x(i) xc],[t_start(i) t]);
            line([x(i-1) xc],[t_start(i-1) t]);
            drawnow;
            if (Periodic&&(i==2)),
               tcoll_left=collision_time([0.0 s(3)], [t t_start(3)], [xleft x(3)],T+1);
               tcoll_right=collision_time([s(nu-2) s(2)], [t_start(nu-2) t],...
                  [x(nu-2) xright],T+1);
               x=[xleft x(3:nu-2) xright xright];
               s=[0 s(3:nu-2) s(2) 0];
               t_start=[t t_start(3:nu-2) t t];
               u=[u(:,3) u(:,3:nu-1) u(:,3) u(:,3)];
               t_coll=[T tcoll_left t_coll(4:nu-2) tcoll_right T T];
            elseif (Periodic&&(i==nu-1)),
               tcoll_right=collision_time([s(nx-2) 0.0],[t_start(nx-2) t],...
                  [x(nx-2) xright],T+1);
               tcoll_left=collision_time([s(nx-1),s(2)],[t t_start(2)],...
                  [xleft, x(2)],T+1);
               x=[xleft xleft x(2:nx-2) xright];
               t_start=[t t t_start(2:nx-2) t];
               s=[0 s(nx-1) s(2:nx-2) 0];
               u=[u(:,nu-2) u(:,nu-2) u(:,2:nu-2) u(:,nu-2)];
               t_coll=[T,T, tcoll_left,t_coll(3:nu-3),tcoll_right,T];
            else
               [ur sr]=feval(Rsolver,u(:,i-1),u(:,i+1),delta);
               nr=size(ur,2);
               ns=length(sr);
               if ns>0,
                  if i>2,
                     t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                        [x(i-2),xc],T+1);
                     if t_coll(i-1)<t-COLL_TOL,
                        disp([-1 t_coll(i-1) t]);
                        disp([s(i-2) sr(1) t_start(i-2) t]);
                        error('not increasing t_coll');
                     end;
                  end;
                  if i<nu-1,
                     t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                        [xc,x(i+1)],T+1);
                     if t_coll(i+1)<=t-COLL_TOL,
                        display([1 t_coll(i+1) t]);
                        display([sr(ns) s(i+1)]);
                        line([xc,xc+sr(ns)*(T-t)],[t,T]);
                        drawnow;
                     end;
                  end;
               else
                  ur=[]; sr=[]; nr=0;
                  if (i>2)&&(i<nu-1),
                     t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                        [t_start(i-2),t_start(i+1)],...
                        [x(i-2),x(i+1)],T+1);
                  end;
               end;
               s=[s(1:i-2) sr s(i+1:nu-1)];
               hone=ones([1,ns]);
               x=[x(1:i-2) xc*hone x(i+1:nu-1)];
               t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
               if ns,
                  u=[u(:,1:i-1) ur u(:,i+1:nu)];
                  t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
               else
                  u=[u(:,1:i-1) u(:,i+2:nu)];
                  t_coll=[t_coll(1:i-1) t_coll(i+2:nu)];
               end;
            end;
            [t i]=min(t_coll);
         end;
         n=length(x);
         for i=1:n
            line([x(i),x(i)+s(i)*(T-t_start(i))],[t_start(i) T]);
         end;
         drawnow;
         x=x+s.*(T-t_start);
         hold off;
      else
         %% Not saving the graphics
         if ~Periodic,
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime);
         else
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime,'periodic',xleft,xright);
         end;
         t_start=startTime*ones(size(x));
         %[xlower xhigher]=get(plotAxis,'XLim');
         first=1;
         [t i]=min(t_coll);
         ncollision=0;
         while (t<T)&&~isempty(s),
            % 				if (~check_states(x,s,t_start,t))
            % 					display('not increasing collisions');
            % 					ncollision
            %				end;
            ncollision=ncollision+1;
            nu=size(u,2);
            nx=length(x);
            if nx~=nu-1,
               display([nx nu]);
               error('Wrong length of x vs u in FrontTrack');
            end;
            xc=x(i)+s(i)*(t-t_start(i));
            % 				if (i>1),
            % 					xc2=x(i-1)+(t-t_start(i-1))*s(i-1);
            % 					if (abs(xc-xc2)>1e-6)
            % 						display(' not same collision point...');
            % 					end;
            % 				end;
            if first,
               first=0;
               h=line([x(i) xc],[t_start(i) t],'EraseMode','none');
            else
               set(h,'XData',[x(i) xc],'YData',[t_start(i) t]);
            end;
            drawnow;
            if i>1,
               set(h,'XData',[x(i-1) xc],'YData',[t_start(i-1) t]);
            end;
            drawnow;
            if (Periodic&&(i==2)),
               tcoll_left=collision_time([0.0 s(3)], [t t_start(3)], [xleft x(3)],T+1);
               tcoll_right=collision_time([s(nu-2) s(2)], [t_start(nu-2) t],...
                  [x(nu-2) xright],T+1);
               x=[xleft x(3:nu-2) xright xright];
               s=[0 s(3:nu-2) s(2) 0];
               t_start=[t t_start(3:nu-2) t t];
               u=[u(:,3) u(:,3:nu-1) u(:,3) u(:,3)];
               t_coll=[T+1 tcoll_left t_coll(4:nu-2) tcoll_right T+1 T+1];
            elseif (Periodic&&(i==nu-1)),
               tcoll_right=collision_time([s(nx-2) 0.0],[t_start(nx-2) t],...
                  [x(nx-2) xright],T+1);
               tcoll_left=collision_time([s(nx-1),s(2)],[t t_start(2)],...
                  [xleft, x(2)],T+1);
               x=[xleft xleft x(2:nx-2) xright];
               t_start=[t t t_start(2:nx-2) t];
               s=[0 s(nx-1) s(2:nx-2) 0];
               u=[u(:,nu-2) u(:,nu-2) u(:,2:nu-2) u(:,nu-2)];
               t_coll=[T+1,T+1, tcoll_left,t_coll(3:nu-3),tcoll_right,T+1];
            else
               [ur sr]=feval(Rsolver,u(:,i-1),u(:,i+1),delta);
               nr =size(ur,2);
               ns=length(sr);
               if ns>0,
                  if i>2,
                     t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                        [x(i-2),xc],T+1);
                     if t_coll(i-1)<t-COLL_TOL,
                        display([-1 t_coll(i-1) t t-t_coll(i-1)]);
                        display([s(i-2) sr(1) t_start(i-2) t]);
                        error('not increasing t_coll');
                     end;
                  end;
                  if i<nu-1,
                     t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                        [xc,x(i+1)],T+1);
                     if t_coll(i+1)<t-COLL_TOL,
                        display([1 t_coll(i+1) t t-t_coll(i+1)]);
                        display([sr(ns) s(i+1)]);
                        error('not increasing t_coll');
                     end;
                  end;
               else
                  ur=[]; sr=[]; nr=0;
                  if (i>2)&&(i<nu-1),
                     t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                        [t_start(i-2),t_start(i+1)],...
                        [x(i-2),x(i+1)],T+1);
                  end;
               end;
               s=[s(1:i-2) sr s(i+1:nu-1)];
               hone=ones([1,ns]);
               x=[x(1:i-2) xc*hone x(i+1:nu-1)];
               t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
               if ns,
                  u=[u(:,1:i-1) ur u(:,i+1:nu)];
                  t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
               else
                  u=[u(:,1:i-1) u(:,i+2:nu)];
                  t_coll=[t_coll(1:i-1) t_coll(i+2:nu)];
               end;
            end;
            [t i]=min(t_coll);
         end;
         n=length(x);
         if first,
            first=0;
            h=line((x(1)+s(1)*(T-t_start(1))),[t_start(1) T],'EraseMode','none');
            drawnow;
         end;
         for i=1:n
            set(h,'XData',[x(i),x(i)+s(i)*(T-t_start(i))],'YData',[t_start(i) T]);
            drawnow;
         end;
         x=x+s.*(T-t_start);
         hold off;
      end;
   elseif polygons,
      %% Using polygons to plot
      if savegraphics,
         
         % ............For drawing filled polygon ...............
         if cmap_set,
            j=cmap;
         else
            j=jet(256);
         end;
         colormap(j);
         NMAP=length(j);
         % ......................................................
         if ~Periodic,
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime);
         else
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime,'periodic',xleft,xright);
         end;
         if (~fixedcolinterval)
            cf=feval(colfunc,u);
            pmax=max(cf);pmin=min(cf);
            clear('cf');
         end;
         t_start=startTime*ones(size(x));
         lims=get(plotAxis,'XLim');
         xlower=lims(1); xhigher=lims(2);
         hold on;
         [t i]=min(t_coll);
         first=1;
         while (t<T)&&~isempty(s),
            nu=size(u,2);
            nx=length(x);
            if nx~=nu-1,
               display([nx nu]);
               error('Wrong length of x vs u in FrontTrack');
            end;
            % ............For drawing filled polygon ...............
            xc=x(i)+s(i)*(t-t_start(i));
            xp=[x(i) xc x(i-1)];
            tp=[t_start(i) t t_start(i-1)];
            colind=clindex(u(:,i),pmax,pmin,NMAP,colfunc);
            patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
            drawnow;
            if i>2,
               xp=[x(i-2) x(i-1) xc];
               tp=[t_start(i-2) t_start(i-1) t];
            else
               xp=[xlower x(i-1) xc xlower];
               tp=[t_start(i-1) t_start(i-1) t t];
            end;
            colind=clindex(u(:,i-1),pmax,pmin,NMAP,colfunc);
            patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
            if i<nu-2,
               xp=[x(i) x(i+1) xc];
               tp=[t_start(i) t_start(i+1) t];
            else
               xp=[x(i) xhigher xhigher xc];
               tp=[t_start(i) t_start(i) t t];
            end;
            colind=clindex(u(:,i+1),pmax,pmin,NMAP,colfunc);
            patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
            % ......................................................
            drawnow;
            if (Periodic&&(i==2)),
               xp=[x(nx-1) x(nx) x(nx)];
               tp=[t_start(nx-1) t_start(nx) t];
               colind=clindex(u(:,nx),pmax,pmin,NMAP,colfunc);
               patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
               drawnow;
               tcoll_left=collision_time([0.0 s(3)], [t t_start(3)], [xleft x(3)],T+1);
               tcoll_right=collision_time([s(nu-2) s(2)], [t_start(nu-2) t],...
                  [x(nu-2) xright],T+1);
               x=[xleft x(3:nu-2) xright xright];
               s=[0 s(3:nu-2) s(2) 0];
               t_start=[t t_start(3:nu-2) t t];
               u=[u(:,3) u(:,3:nu-1) u(:,3) u(:,3)];
               t_coll=[T tcoll_left t_coll(4:nu-2) tcoll_right T+1 T+1];
            elseif (Periodic&&(i==nu-1)),
               xp=[x(1) x(1) x(2)];
               tp=[t_start(1) t t_start(2)];
               colind=clindex(u(:,2),pmax,pmin,NMAP,colfunc);
               patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
               tcoll_right=collision_time([s(nx-2) 0.0],[t_start(nx-2) t],...
                  [x(nx-2) xright],T+1);
               tcoll_left=collision_time([s(nx-1),s(2)],[t t_start(2)],...
                  [xleft, x(2)],T+1);
               x=[xleft xleft x(2:nx-2) xright];
               t_start=[t t t_start(2:nx-2) t];
               s=[0 s(nx-1) s(2:nx-2) 0];
               u=[u(:,nu-2) u(:,nu-2) u(:,2:nu-2) u(:,nu-2)];
               t_coll=[T+1,T+1, tcoll_left,t_coll(3:nu-3),tcoll_right,T+1];
            else
               [ur sr]=feval(Rsolver,u(:,i-1),u(:,i+1),delta);
               nr=size(ur,2);
               ns=length(sr);
               if ns>0,
                  if i>2,
                     t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                        [x(i-2),xc],T+1);
                     if t_coll(i-1)<t-COLL_TOL,
                        display([-1 t_coll(i-1) t t-t_coll(i-1)]);
                        display([s(i-2) sr(1) t_start(i-2) t]);
                        error('not increasing t_coll');
                     end;
                  end;
                  if i<nu-1,
                     t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                        [xc,x(i+1)],T+1);
                     if t_coll(i+1)<t-COLL_TOL,
                        display([1 t_coll(i+1) t t-t_coll(i+1)]);
                        display([sr(ns) s(i+1)]);
                        error('not increasing t_coll');
                     end;
                  end;
               else
                  ur=[]; sr=[]; nr=0;
                  if (i>2)&&(i<nu-1),
                     t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                        [t_start(i-2),t_start(i+1)],...
                        [x(i-2),x(i+1)],T+1);
                  end;
               end;
               s=[s(1:i-2) sr s(i+1:nu-1)];
               hone=ones([1,ns]);
               x=[x(1:i-2) xc*hone x(i+1:nu-1)];
               t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
               if ns,
                  u=[u(:,1:i-1) ur u(:,i+1:nu)];
                  t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
               else
                  u=[u(:,1:i-1) u(:,i+2:nu)];
                  t_coll=[t_coll(1:i-1) t_coll(i+2:nu)];
               end;
            end;
            [t i]=min(t_coll);
         end;
         %     For drawing filled polygon
         xp=[xlower xlower];
         tp=[t_start(1) T];
         for i=1:n
            xp2=[x(i)+s(i)*(T-t_start(i)) x(i)];
            tp2=[T t_start(i)];
            colind=clindex(u(:,i),pmax,pmin,NMAP,colfunc);
            patch([xp xp2],[tp tp2],j(colind,:),'EdgeColor','none','EraseMode','none');
            xp=xp2(2:-1:1); tp=tp2(2:-1:1);
         end;
         drawnow;
         i=n;
         xp2=[xhigher xhigher]; tp2=[T t_start(i)];
         colind=clindex(u(:,i+1),pmax,pmin,NMAP,colfunc);
         patch([xp xp2],[tp tp2],j(colind,:),'EdgeColor','none','EraseMode','none');
         %
         x=x+s.*(T-t_start);
         hold off;
      else
         %% Polygons without graphics save
         %       For drawing filled polygon
         if cmap_set,
            j=cmap;
         else
            j=jet(256);
         end;
         colormap(j);
         NMAP=length(j);
         %
         if ~Periodic,
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime);
         else
            [u s t_coll,x]=initialize(u0,x0,T,delta,'Riemann_solver',Rsolver,...
               'start_time',startTime,'periodic',xleft,xright);
         end;
         if (~fixedcolinterval)
            cf=feval(colfunc,u);
            pmax=max(cf);pmin=min(cf);
            clear('cf');
         end;
         t_start=startTime*ones(size(x));
         lims=get(plotAxis,'XLim');
         xlower=lims(1); xhigher=lims(2);
         hold on;
         [t i]=min(t_coll);
         first=1;
         while (t<T)&&~isempty(s),
            nu=size(u,2);
            nx=length(x);
            if nx~=nu-1,
               display([nx nu]);
               error('Wrong length of x vs u in FrontTrack');
            end;
            %  For drawing filled polygon
            xc=x(i)+s(i)*(t-t_start(i));
            xp=[x(i) xc x(i-1)];
            tp=[t_start(i) t t_start(i-1)];
            colind=clindex(u(:,i),pmax,pmin,NMAP,colfunc);
            if first,
               h=patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
               first=0;
            else
               set(h,'XData',xp,'YData',tp,'FaceColor',j(colind,:));
            end;
            drawnow;
            if i>2,
               xp=[x(i-2) x(i-1) xc];
               tp=[t_start(i-2) t_start(i-1) t];
            else
               xp=[xlower x(i-1) xc xlower];
               tp=[t_start(i-1) t_start(i-1) t t];
            end;
            colind=clindex(u(:,i-1),pmax,pmin,NMAP,colfunc);
            set(h,'XData',xp,'YData',tp,'FaceColor',j(colind,:));
            %	drawnow;
            if i<nu-2,
               xp=[x(i) x(i+1) xc];
               tp=[t_start(i) t_start(i+1) t];
            else
               xp=[x(i) xhigher xhigher xc];
               tp=[t_start(i) t_start(i) t t];
            end;
            colind=clindex(u(:,i+1),pmax,pmin,NMAP,colfunc);
            set(h,'XData',xp,'YData',tp,'FaceColor',j(colind,:));
            drawnow;
            if (Periodic&&(i==2)),
               xp=[x(nx-1) x(nx) x(nx)];
               tp=[t_start(nx-1) t_start(nx) t];
               colind=clindex(u(:,nx),pmax,pmin,NMAP,colfunc);
               if first,
                  first=0;
                  h=patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
               else
                  set(h,'XData',xp,'YData',tp,'FaceColor',j(colind,:));
               end;
               tcoll_left=collision_time([0.0 s(3)], [t t_start(3)], [xleft x(3)],T+1);
               tcoll_right=collision_time([s(nu-2) s(2)], [t_start(nu-2) t],...
                  [x(nu-2) xright],T+1);
               x=[xleft x(3:nu-2) xright xright];
               s=[0 s(3:nu-2) s(2) 0];
               t_start=[t t_start(3:nu-2) t t];
               u=[u(:,3) u(:,3:nu-1) u(:,3) u(:,3)];
               t_coll=[T tcoll_left t_coll(4:nu-2) tcoll_right T T];
            elseif (Periodic&&(i==nu-1)),
               xp=[x(1) x(1) x(2)];
               tp=[t_start(1) t t_start(2)];
               colind=clindex(u(:,2),pmax,pmin,NMAP,colfunc);
               if first,
                  first=0;
                  h=patch(xp,tp,j(colind,:),'EdgeColor','none','EraseMode','none');
               else
                  set(h,'XData',xp,'YData',tp,'FaceColor',j(colind,:));
               end;
               tcoll_right=collision_time([s(nx-2) 0.0],[t_start(nx-2) t],...
                  [x(nx-2) xright],T+1);
               tcoll_left=collision_time([s(nx-1),s(2)],[t t_start(2)],...
                  [xleft, x(2)],T+1);
               x=[xleft xleft x(2:nx-2) xright];
               t_start=[t t t_start(2:nx-2) t];
               s=[0 s(nx-1) s(2:nx-2) 0];
               u=[u(:,nu-2) u(:,nu-2) u(:,2:nu-2) u(:,nu-2)];
               t_coll=[T,T, tcoll_left,t_coll(3:nu-3),tcoll_right,T];
            else
               [ur sr]=feval(Rsolver,u(:,i-1),u(:,i+1),delta);
               nr=size(ur,2);
               ns=length(sr);
               if ns>0,
                  if i>2,
                     t_coll(i-1)=collision_time([s(i-2) sr(1)],[t_start(i-2),t],...
                        [x(i-2),xc],T+1);
                     if t_coll(i-1)<t-COLL_TOL,
                        display([-1 t_coll(i-1) t t-t_coll(i-1)]);
                        display([s(i-2) sr(1) t_start(i-2) t]);
                        error('not increasing t_coll');
                     end;
                  end;
                  if i<nu-1,
                     t_coll(i+1)=collision_time([sr(ns),s(i+1)],[t,t_start(i+1)],...
                        [xc,x(i+1)],T+1);
                     if t_coll(i+1)<t-COLL_TOL,
                        display([1 t_coll(i+1) t t-t_coll(i+1)]);
                        display([sr(ns) s(i+1)]);
                        error('not increasing t_coll');
                     end;
                  end;
               else
                  ur=[]; sr=[]; nr=0;
                  if (i>2)&&(i<nu-1),
                     t_coll(i-1)=collision_time([s(i-2) s(i+1)],...
                        [t_start(i-2),t_start(i+1)],...
                        [x(i-2),x(i+1)],T+1);
                  end;
               end;
               s=[s(1:i-2) sr s(i+1:nu-1)];
               hone=ones([1,ns]);
               x=[x(1:i-2) xc*hone x(i+1:nu-1)];
               t_start=[t_start(1:i-2) t*hone t_start(i+1:nu-1)];
               if ns,
                  u=[u(:,1:i-1) ur u(:,i+1:nu)];
                  t_coll=[t_coll(1:i-1) (T+1)*ones([1,nr]) t_coll(i+1:nu)];
               else
                  u=[u(:,1:i-1) u(:,i+2:nu)];
                  t_coll=[t_coll(1:i-1) t_coll(i+2:nu)];
               end;
            end;
            [t i]=min(t_coll);
         end;
         % For drawing filled polygon
         n=length(x);
         xp=[xlower xlower];
         tp=[t_start(1) T];
         for i=1:n
            xp2=[x(i)+s(i)*(T-t_start(i)) x(i)];
            tp2=[T t_start(i)];
            colind=clindex(u(:,i),pmax,pmin,NMAP,colfunc);
            if first,
               first=0;
               h=patch([xp xp2],[tp tp2],'Facecolor',j(colind,:),'EraseMode','none');
            else
               set(h,'Xdata',[xp xp2],'YData',[tp tp2],'FaceColor',j(colind,:));
            end;
            xp=xp2(2:-1:1); tp=tp2(2:-1:1);
         end;
         drawnow;
         i=n;
         xp2=[xhigher xhigher]; tp2=[T t_start(i)];
         colind=clindex(u(:,i+1),pmax,pmin,NMAP,colfunc);
         set(h,'Xdata',[xp xp2],'YData',[tp tp2],'FaceColor',j(colind,:));
         drawnow;
         x=x+s.*(T-t_start);
         hold off;
      end;
   end;
end;

%% Apology
%  This code was started long before  matlab had built in objects.
%  Therefore the code is based on arrays which constantly have to be
%  changed. This takes a long time, and using pointers would probably
%  increase the speed of the code. For "serious" computations I use a
%  C-code with pointers.  -nilshr
%
