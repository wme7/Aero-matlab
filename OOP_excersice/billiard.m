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
            %whatever.XMIN = a;
            %whatever.XMAX = b;
            %whatever.YMIN = c;
            %whatever.YMAX = d;
            
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
          
          % draw a single ball
          x = 10;
          y = 25;
          
          r = rectangle('position',[x,y,w/20,w/20],...
                    'FaceColor','b',...
                    'curvature',1.0 ...
                    );
           
           % the initial position
           x_new = x;
           y_new = y;
           
           % update position
           
           for time = 1:450;
           
           x_new = x_new+0.15;
           y_new = 20*sin(4*x_new/h);
           
           r.Position = [x_new,y_new+25,w/20,w/20];
           
           drawnow
           
           end
           
        end
          
    end
end