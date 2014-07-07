function u0 = IC2d(x,y,option)
%% Initial Condition (IC)
% This subroutine creates an special IC for testing purposes.
% The present IC models a rectangular domain with four equally sized
% regions with diferrent initial state value, u, an adimensional value.
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06

%% Select your IC

switch option
    case 1 % 4 quadrants in Domain [-1,1]^2
        % Preallocate u0,
        u0 = zeros(size(x));

        % Initial Condition for our 2D domain
        u0(x >0 & y >0) = 0.50; % region 1
        u0(x<=0 & y >0) = 0.70; % region 2
        u0(x<=0 & y<=0) = 0.10; % region 3
        u0(x >0 & y<=0) = 0.90; % region 4
        
    case 2 % Square Jump
        % set center
        x0 = 0.0; y0 = 0.0;
        
        % parameters
        Dx = x(end)-x(1);
        Dy = y(end)-y(1);
        
        % Preallocate u0
        u0 = ones(size(x));
        
        % Parameters of region
        x1 = x0+Dx/10; x2 = x0-Dx/10;
        y1 = y0+Dy/10; y2 = y0-Dy/10;
        
        % Build Jump
        u0(x>x1) = 0.1;     u0(x<x2) = 0.1;
        u0(y>y1) = 0.1;     u0(y<y2) = 0.1;
        
    case 3 % Sine*Cosine 2-D in Domain [-1,1]^2
        u0 = sin(pi*x).*cos(pi*y);
        
    case 4 % Gaussian Jump
        % set center
        x0 = 0.0; y0 = 0.0; 
        
        % Gaussian
        u0 = 0.1 + 0.5*exp(-20*((x-x0).^2+(y-y0).^2));
        
    case 5 % Cylindrical Jump
        % set center
        x0 = 0.0; y0 = 0.0;
        
        % radious
        r = 0.5;
        
        % Gaussian
        u0 = 0.1*ones(size(x));
        u0(sqrt((x+x0).^2+(y+y0).^2)<r) = 1.0;
        
    case 6 % rectangle in x direction
        % set center
        x0 = 0.0;
        
        % parameters
        Dx = x(end)-x(1);
                
        % Preallocate u0
        u0 = ones(size(x));
        
        % Parameters of region
        x1 = x0+Dx/12; x2 = x0-Dx/12;
                
        % Build Jump
        u0(x>x1) = 0.1;     u0(x<x2) = 0.1;
                
	case 7 % rectangle in y direction
        % set center
        y0 = 0.0;
        
        % parameters
        Dy = y(end)-y(1);
        
        % Preallocate u0
        u0 = ones(size(x));
        
        % Parameters of region
        y1 = y0+Dy/12; y2 = y0-Dy/12;
        
        % Build Jump
        u0(y>y1) = 0.1;     u0(y<y2) = 0.1;
    case 8 % Riemann for range [-1,1]
        x1=x<0.0; x2=x>=0.0; 
        
        % Riemann IC
        u0=1*x1+0.1*x2;
        
    otherwise 
        error('case not listed :P')
end
