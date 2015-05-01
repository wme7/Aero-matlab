classdef billiard
    
    % this are the dimensions of the table
    properties(Constant = true)
        XMAX = 100;
        XMIN = 0;
        YMAX = 50;
        YMIN = 0;
    end

    methods
        function whatever = billiard()
            % step 1. Draw the billiard table
            whatever.drawtable();
         
            % step 2. Draw the balls
          % ===========================
          % draw a single ball
          x = 10;   % Initial x-position
          y = 25;   % Initial y-position
          u = 1;    % Initial x-velocity (moves in x # pixels per iteration)
          v = 1;    % Initial y-velcotiy (moves in y # pixels per iteration)
          d = 5;    % Diameter of the ball
          m = 1;    % The mass of the ball
          
          ball = rectangle('position',[x,y,d,d],...
                    'FaceColor','b',...
                    'curvature',1.0 ...
                    );
          % ===========================
            % Update position
          % ===========================
           % the initial position
           x_new = x;
           y_new = y;
           
           % Update position
           
           dt = 1; % the time step is constant
           
           for iteration = 1:950;
               
               disp(x_new)
               
               if (abs(x_new-whatever.XMAX+d)<eps) || (abs(x_new-whatever.XMIN)<eps)
                   u = -u;
               end
               
               if (abs(y_new-whatever.YMAX+d)<eps) || (abs(y_new-whatever.YMIN)<eps)
                   v = -v;
               end
               
           x_new = x_new + u*dt;
           y_new = y_new + v*dt;
           
           ball.Position = [x_new,y_new,d,d];
           
           drawnow
           
           end

           % ===========================
        end
        
        function drawtable(whatever)
            % properties of the table
            w = whatever.XMAX-whatever.XMIN;            
            h = whatever.YMAX-whatever.YMIN;
            x0 = whatever.XMIN;
            y0 = whatever.YMIN;
            % draw the table
            figure( 'name'   ,'TheCrazyBilliard2000',...
                    'menubar','none',...
                    'units'  ,'pixels'...
                    );
            rectangle('position',[x0,y0,w,h],...
                      'FaceColor','w',...
                      'curvature',0.1...
                      );
            axis('equal','off');
        end
          
    end
end