classdef Balls<handle

%--------------------------------------------------------------------------
%Collision Ball Libraries
%Version 1.00
%Created by Stepen
%Created 20 March 2013
%Last modified 20 May 2013
%--------------------------------------------------------------------------
%Constructor Format:
%- Balls()
%  creates an empty Collision Ball object.
%- Balls(x,y,u,v)
%  creates a Collision Ball where x, y, u, and v respectively are column
%  array with same row length defining balls position and velocity.
%- Balls(x,y,u,v,m,r,c)
%  creates a Collision Ball where x, y, u, v, m, r, and c respectively are
%  column array (except for c is 3column array) with same row length
%  defining balls position, velocity, mass, radius, and color.
%Public Non-static Method:
%- addBall(x,y,u,v,m,r,c)
%  add a ball to Collision Ball object.
%- addRandomBall(x,y,u,v,m,r,c)
%  add a random ball to Collision Ball object.
%- removeBall(ballID)
%  remove selected ball from Collision Ball object.
%- isOccupied(x,y,r)
%  check whether a ball at x and y position with radius of r can be added
%  to Collision Ball ojbect without conflict.
%- moveBall(dt)
%  move all balls in Collision Ball object for given time interval.
%- draw()
%  plot all balls in Collision Ball object at current time.
%- play(dt)
%  move and animate all balls in Collision Ball object using given time
%  interval as animation refresh rate.
%Public Static Method:
%- isOutBound(x,y,radius)
%  check a ball at x and y position with radius of r is within Collision
%  Ball object's containment space.
%--------------------------------------------------------------------------

%CodeStart-----------------------------------------------------------------
%Declaring object properties
    properties(SetAccess=protected)
        n_ball
        mass
        radius
        color
        x
        y
        u
        v
        ANIMATIONSTAT=false;
    end
%Declaring constant
    properties(Constant=true)
        MASS_MAX=100;
        RADIUS_MAX=20;
        XMIN=-100;
        XMAX=100;
        YMIN=-100;
        YMAX=100;
        UMIN=-25;
        UMAX=25;
        VMIN=-25;
        VMAX=25;
    end
%Declaring constructor
    methods
        function this=Balls(varargin)
            %Initiating field
            this.n_ball=0;
            this.mass=zeros(0,1);
            this.radius=zeros(0,1);
            this.color=zeros(0,3);
            this.x=zeros(0,1);
            this.y=zeros(0,1);
            this.u=zeros(0,1);
            this.v=zeros(0,1);
            %Assigning input arguments to its appropriate field
            switch nargin
                case 7
                    for i=1:numel(varargin{1})
                        this.addBall(varargin{1}(i),varargin{2}(i),...
                                     varargin{3}(i),varargin{4}(i),...
                                     varargin{5}(i),varargin{6}(i),...
                                     varargin{7}(i,:))
                    end
                case 4
                    for i=1:numel(varargin{1})
                        this.addBall(varargin{1}(i),varargin{2}(i),...
                                     varargin{3}(i),varargin{4}(i),...
                                     10,10,[0,0,0])
                    end
                case 1
                    Balls.checksinglevarargin(varargin{1});
                    for i=1:varargin{1}
                        ballValidity=false;
                        while ~ballValidity
                            %Creating random ball
                            [xr,yr,ur,vr,mr,rr,cr]=Balls.generateball();
                            %Checking random ball validity
                            ballValidity=(~this.isOccupied(xr,yr,rr))&&...
                                         (~Balls.isOutBound(xr,yr,rr));
                            %Adding ball
                            if ballValidity
                                this.addBall(xr,yr,ur,vr,mr,rr,cr)
                            end
                        end
                    end
                case 0
                    
                otherwise
                    error('Unexpected constructor format!');
            end
        end
    end
%Declaring public method
    methods
        %Declaring method to add a ball
        function addBall(this,x,y,u,v,mass,radius,color)
            %Checking input arguments
            switch nargin
                case 1
                    this.addRandomBall()
                    return;
                case 5
                    Balls.checkfourvarargins(x,u,v,y)
                    mass=10;
                    radius=10;
                    color=[0,0,255];
                case 8
                    Balls.checksevenvarargins(x,y,u,v,mass,radius,color);
                otherwise
                    error('Unexpected format of input argument!');
            end
            %Checking for validity
            if this.isOccupied(x,y,radius)
                error('New ball is conflicting with the existing balls!');
            end
            if Balls.isOutBound(x,y,radius)
                error('New ball is out of bounds!');
            end
            %Adding ball
            this.n_ball=this.n_ball+1;
            this.mass=[this.mass;mass];
            this.radius=[this.radius;radius];
            this.color=[this.color;color];
            this.x=[this.x;x];
            this.y=[this.y;y];
            this.u=[this.u;u];
            this.v=[this.v;v];
            %Updating drawing plot (if drawing already exist)
            mainfig=Balls.getfigurehandle();
            if ~isempty(mainfig)
                rectangle('Parent',Balls.getaxeshandle(),...
                          'FaceColor',color/255,...
                          'Curvature',[1,1],...
                          'Position',[x-radius,y-radius,...
                                      2*radius,2*radius]);
            end
        end
        %Declaring method to add a random ball
        function addRandomBall(this)
            %Initiating variable for while loop
            canBeAdded=false;
            %Generating random ball until it is valid
            while ~canBeAdded
                %Generating random position
                [xr,yr,ur,vr,mr,rr,cr]=Balls.generateball();
                %Checking random ball validity
                canBeAdded=(~this.isOccupied(xr,yr,rr))&&...
                           (~Balls.isOutBound(xr,yr,rr));
            end
            %Adding random ball
            this.addBall(xr,yr,ur,vr,mr,rr,cr);
        end
        %Declaring method to remove a ball
        function removeBall(this,ballID)
            %Checking input value
            Balls.checksinglevarargin(ballID);
            if ballID>this.n_ball
                error('Ball ID should be less that number of balls!')
            end
            %Deleting ball
            this.mass(ballID)=[];
            this.radius(ballID)=[];
            this.color(ballID,:)=[];
            this.x(ballID)=[];
            this.y(ballID)=[];
            this.u(ballID)=[];
            this.v(ballID)=[];
            this.n_ball=this.n_ball-1;
            %Updating ball drawing (if drawing already exist)
            mainfig=Balls.getfigurehandle();
            if ~isempty(mainfig)
                rectlist=findobj(mainfig,'Type','rect');
                delete(rectlist(this.n_ball+2-ballID))
            end
        end
        %Declaring method to check whether a position is available
        function stat=isOccupied(this,x,y,radius)
            %Initiating output variable
            stat=false;
            %Checking space for each ball
            for i=1:this.n_ball
                S=norm([this.x(i)-x,this.y(i)-y]);
                if S<(this.radius(i)+radius)
                    stat=true;
                    return
                end
            end
        end
        %Declaring function to move the balls
        function moveBall(this,dt)
            %Predicting ball position without considering collision
            [xf,yf]=predictposition(this,dt);
            %Finding collision
            stat=Balls.findcollision(xf,yf,this.radius);
            %while ~isempty(stat)
                %Correcting the speed of colliding balls
                for i=1:size(stat,1)
                    ID1=stat(i,1);
                    ID2=stat(i,2);
                    [this.u(ID1),this.v(ID1),...
                     this.u(ID2),this.v(ID2)]=Balls.collide(this.mass(ID1),...
                                                            this.x(ID1),...
                                                            this.y(ID1),...
                                                            this.u(ID1),...
                                                            this.v(ID1),...
                                                            this.mass(ID2),...
                                                            this.x(ID2),...
                                                            this.y(ID2),...
                                                            this.u(ID2),...
                                                            this.v(ID2));
                end
                %Finding out-of-bount balls
                hstat=(xf<Balls.XMIN+this.radius)|...
                      (xf>Balls.XMAX-this.radius);
                this.u(hstat)=-this.u(hstat);
                vstat=(yf<Balls.YMIN+this.radius)|...
                      (yf>Balls.YMAX-this.radius);
                this.v(vstat)=-this.v(vstat);
                %Limiting maximum speed
                this.u(this.u<Balls.UMIN)=Balls.UMIN;
                this.u(this.u>Balls.UMAX)=Balls.UMAX;
                this.v(this.v<Balls.VMIN)=Balls.VMIN;
                this.v(this.v>Balls.VMAX)=Balls.VMAX;
                %Recalculating ball position with updated velocity
                [xf,yf]=predictposition(this,dt);
                %Rechecking collision
                %stat=Balls.findcollision(xf,yf,this.radius);
                %Updating ball position
                this.x=xf;
                this.y=yf;
            %end
        end
        %Declaring function to draw the balls
        function draw(this)
            %Checking for existing figure
            mainfig=Balls.getfigurehandle();
            %Creating or getting figure handle
            if isempty(mainfig)
                %Creating brand new figure and axes
                WinSize=[1,1,800,800];
                ScrSize=get(0,'ScreenSize');
                mainfig=figure('Name','Elastic Collision Balls',...
                               'NumberTitle','off',...
                               'Menubar','none',...
                               'Units','pixels',...
                               'Position',[(ScrSize(3)-WinSize(3))/2,...
                                           (ScrSize(4)-WinSize(4))/2,...
                                           WinSize(3),WinSize(4)]);
                mainaxs=axes('Parent',mainfig,...
                             'Units','normalized',...
                             'Position',[0.1,0.1,0.8,0.8]);
            else
                %Getting axes handle
                mainaxs=Balls.getaxeshandle();
                %Deleting previous rectangle object
                delete(Balls.getrecthandle());
            end
            %Drawing white rectangle as field
            rectangle('Parent',mainaxs,...
                      'FaceColor','w',...
                      'Position',[Balls.XMIN,Balls.YMIN,...
                                  Balls.XMAX-Balls.XMIN,...
                                  Balls.YMAX-Balls.YMIN]);
            %Drawing balls as curve rectangle
            for i=1:this.n_ball
                rectangle('Parent',mainaxs,...
                          'FaceColor',this.color(i,:)/255,...
                          'Curvature',[1,1],...
                          'Position',[this.x(i)-...
                                      this.radius(i),...
                                      this.y(i)-...
                                      this.radius(i),...
                                      2*this.radius(i),...
                                      2*this.radius(i)]);
            end
            %Correcting axes
            axis(mainaxs,[Balls.XMIN,Balls.XMAX,Balls.YMIN,Balls.YMAX]);
            axis(mainaxs,'equal');
            axis(mainaxs,'off');
        end
        %Declaring function to animate the balls
        function play(this,dt)
            %Checking for existing figure
            mainfig=Balls.getfigurehandle();
            %Creating or getting figure handle
            if isempty(mainfig)
                %Creating brand new figure and axes
                WinSize=[1,1,800,800];
                ScrSize=get(0,'ScreenSize');
                mainfig=figure('Name','Elastic Collision Balls',...
                               'NumberTitle','off',...
                               'Menubar','none',...
                               'Units','pixels',...
                               'Position',[(ScrSize(3)-WinSize(3))/2,...
                                           (ScrSize(4)-WinSize(4))/2,...
                                           WinSize(3),WinSize(4)]);
                mainaxs=axes('Parent',mainfig,...
                             'Units','normalized',...
                             'Position',[0.1,0.1,0.8,0.8]);
            else
                %Getting axes handle
                mainaxs=Balls.getaxeshandle();
                %Deleting previous rectangle object
                delete(Balls.getrecthandle());
            end
            %Drawing rectangle for visualization field
            rectangle('Parent',mainaxs,...
                      'FaceColor','w',...
                      'Position',[Balls.XMIN,Balls.YMIN,...
                                  Balls.XMAX-Balls.XMIN,...
                                  Balls.YMAX-Balls.YMIN]);
            %Drawing balls as curved rectangle
            ballrect=zeros(this.n_ball,1);
            for i=1:this.n_ball
                ballrect(i)=rectangle('Parent',mainaxs,...
                                      'FaceColor',this.color(i,:)/255,...
                                      'Curvature',[1,1],...
                                      'Position',[this.x(i)-...
                                                  this.radius(i),...
                                                  this.y(i)-...
                                                  this.radius(i),...
                                                  2*this.radius(i),...
                                                  2*this.radius(i)]);
            end
            %Animating balls
            this.ANIMATIONSTAT=true;
            while this.ANIMATIONSTAT
                %Timing process for animation delay
                t=tic;
                %Moving balls
                moveBall(this,dt)
                %Redrawing balls
                for i=1:this.n_ball
                    try
                        set(ballrect(i),'Position',[this.x(i)-...
                                                    this.radius(i),...
                                                    this.y(i)-...
                                                    this.radius(i),...
                                                    2*this.radius(i),...
                                                    2*this.radius(i)])
                    catch
                        this.ANIMATIONSTAT=false;
                        return
                    end
                end
                %Correcting axes
                axis(mainaxs,[Balls.XMIN,Balls.XMAX,Balls.YMIN,Balls.YMAX]);
                axis(mainaxs,'equal');
                axis(mainaxs,'off');
                %Delaying animation
                tt=toc(t);
                pause(dt-tt)
                if (dt-tt)<0
                    disp('Interval time too short for animation!')
                    disp('Animation might not be displayed properly!')
                end
            end
        end
    end
%Declaring public static method
    methods(Static=true)
    %Declaring method to check whether a ball is out of bound
        function stat=isOutBound(x,y,radius)
            %Initiating output variable
            stat=false;
            %Checking space for ball position
            if (x<Balls.XMIN+radius)||(x>Balls.XMAX-radius)||...
               (y<Balls.YMIN+radius)||(y>Balls.YMAX-radius)
                stat=true;
                return;
            end
        end
    end
%Declaring private method
    methods(Access=protected)
        %Declaring function to move ball by neglecting collision
        function [x,y]=predictposition(this,dt)
            x=this.x+(this.u*dt)/2;
            y=this.y+(this.v*dt)/2;
        end
    end
%Declaring static private methods
    methods(Access=protected,...
            Static=true)
        %Declaring function to get Balls' figure
        function fig=getfigurehandle()
            %Initiating output variable
            fig=[];
            %Finding all figure objects
            list=findobj(0,'type','figure');
            %Checking each figure name
            for i =1:numel(list)
                if isequal(get(list(i),'Name'),'Elastic Collision Balls')
                    fig=list(i);
                    return;
                end
            end
        end
        %Declaring object to get Balls' axes
        function axs=getaxeshandle()
            axs=findobj(Balls.getfigurehandle(),'type','axes');
        end
        %Declaring function to get Balls' rectangle object
        function rectlist=getrecthandle()
            rectlist=findobj(Balls.getfigurehandle(),'type','rect');
        end
        %Declaring function to check input argument for constructor
        function checksinglevarargin(var1)
            %Checking varargin array size
            if numel(var1)~=1
                error('Input number of balls must be a scalar!');
            end
            %Checking varargin data type
            if (var1<=0)||(mod(var1,1)~=0)
                error('Input number of balls must be a positive integer!');
            end
        end
        function checkfourvarargins(var1,var2,var3,var4)
            %Checking varargin array size
            if (numel(var1)~=1)||(numel(var2)~=1)||...
               (numel(var3)~=1)||(numel(var4)~=1)
                error('Input x, y, u, and v must be a scalar!');
            end
            %Checking varargin data type
            if (~isreal(var1))||(~isreal(var2))||...
               (~isreal(var3))||(~isreal(var4))
                error('Input x, y, u, and v must be real numbers!');
            end
        end
        function checksevenvarargins(var1,var2,var3,var4,var5,var6,var7)
            %Checking varargin for x,y,u,v
            Balls.checkfourvarargins(var1,var2,var3,var4);
            %Checking varargin dimension size
            if (numel(var5)~=1)||(numel(var6)~=1)||(numel(var7)~=3)
                error(['Input mass, radius, and color must be a 2D',...
                       ' column array (3 column array for color)!']);
            end
            %Checking varargin data type
            if (~isreal(var5))||(var5<=0)%||(var5>Balls.MASS_MAX)
                error(['Input mass must be a positive real number',...
                       ' smaller than ',num2str(Balls.MASS_MAX),'!']);
            end
            if (~isreal(var6))||(var6<=0)||(var6>Balls.RADIUS_MAX)
                error(['Input radius must be a positive real number',...
                       ' smaller than ',num2str(Balls.RADIUS_MAX),'!']);
            end
            if (sum(var7<0)~=0)||(sum(var7>255)~=0)||...
               (sum(mod(var7,1)~=0)~=0)
                error(['Input color must be a positive integer',...
                       ' between 0 and 255!']);
            end
        end
        %Declaring function to generate random ball
        function [xr,yr,ur,vr,mr,rr,cr]=generateball()
            mr=rand*Balls.MASS_MAX;
            rr=rand*Balls.RADIUS_MAX;
            cr=floor(rand(1,3)*256);
            xr=Balls.XMIN+rr+...
                (rand*(Balls.XMAX-Balls.XMIN-(2*rr)));
            yr=Balls.YMIN+rr+...
                (rand*(Balls.YMAX-Balls.XMIN)-(2*rr));
            ur=Balls.UMIN+(rand*(Balls.UMAX-Balls.UMIN));
            vr=Balls.VMIN+(rand*(Balls.VMAX-Balls.VMIN));
        end
        %Declaring function to find colliding balls
        function stat=findcollision(x,y,r)
            %Preallocating array for output variable
            stat=zeros(0,2);
            %Finding collision
            for i=1:numel(x)-1
                if isempty(find(stat(:,2)==i, 1))
                    for j=i+1:numel(x)
                        S=norm([x(i)-x(j),y(i)-y(j)]);
                        if S<(r(i)+r(j))
                            stat=[stat;[i,j]];
                            %break
                        end
                    end
                end
            end
        end
        %Declaring function to resolve 2D collision
        function [uf1,vf1,uf2,vf2]=collide(m1,x1,y1,u1,v1,m2,x2,y2,u2,v2)
            %Finding normal angle between two balls
            theta=atan2(y2-y1,x2-x1);
            %Transforming balls velocity to normal coordinate
            VN1=(u1*cos(theta))+(v1*sin(theta));
            VT1=(-u1*sin(theta))+(v1*cos(theta));
            VN2=(u2*cos(theta))+(v2*sin(theta));
            VT2=(-u2*sin(theta))+(v2*cos(theta));
            %Resolving colllision for normal axis
            VN1_F=((VN1*(m1-m2))+(2*m2*VN2))/(m1+m2);
            VN2_F=((VN2*(m2-m1))+(2*m1*VN1))/(m1+m2);
            %Retransforming balls velocity to original coordinate
            uf1=(VN1_F*cos(theta))-(VT1*sin(theta));
            vf1=(VN1_F*sin(theta))+(VT1*cos(theta));
            uf2=(VN2_F*cos(theta))-(VT2*sin(theta));
            vf2=(VN2_F*sin(theta))+(VT2*cos(theta));
        end
    end
%CodeEnd-------------------------------------------------------------------
    
end