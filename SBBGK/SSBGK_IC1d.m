function [r_0,u_0,t_0,p_0,rho_0,E_0,tEnd,cfl] = SSBGK_IC1d(x,input)
% Load the IC of a classical 2D Riemann Problem configuration. 
% By Manuel Diaz 2012.10.24.
% In the notation we take advantage of the matlab array notation as follows
%
% prop = [prop_left , prop_right]
% 
% Notation:
% u = Velocity in x direction
% p = Pressure
% rho = Density
% r = Fugacity
% E = Energy
% t = temperature
%
% Based on:
% http://wonka.physics.ncsu.edu/pub/VH-1/testpage/ and
% http://sitemaker.umich.edu/anand/files/riemann_shock-tube.pdf
% See also my routine CFD/Riemann.m to compute exact solutions.
%
%% Initial Physical Properties per case:
switch input
    case{1} % Configuration 1, Sod's Problem
        fprintf('Case 1: Sods problem \n');
        p   = [1    0.1  ];
        u   = [0.75 0    ];
        rho = [1    0.125];
        tEnd = 0.1; cfl = 0.15;
        
    case{2} % Configuration 2, Left Expansion and right strong shock
        fprintf('Case 2: Left Expansion and right strong shock \n');
        p   = [1000 0.1  ];
        u   = [0    0    ];
        rho = [3    2    ];
        tEnd = 0.02; cfl = 0.15;
        
    case{3} % Configuration 3, Right Expansion and left strong shock
        fprintf('Case 3: Right Expansion and left strong shock \n');
        p   = [7    10   ];
        u   = [0    0    ];
        rho = [1    1    ];
        tEnd = 0.1; cfl = 0.15;
        
    case{4} % Configuration 4, Double Shock
        fprintf('Case 4: Double Shock \n');
        p   = [450  45   ];
        u   = [20   -6   ];
        rho = [6    6    ];
        tEnd = 0.02; cfl = 0.15;
        
    case{5} % Configuration 5, Double Expansion
        fprintf('Case 5: Double Expansion \n');
        p   = [40   40   ];
        u   = [-2   2    ];
        rho = [1    2.5  ];
        tEnd = 0.03; cfl = 0.15;

    case{6} % Configuration 6, Cavitation
        fprintf('Case 6: Cavitation \n');
        p   = [0.4  0.4  ];
        u   = [-20  20   ];
        rho = [1    1    ];
        tEnd = 0.1; cfl = 0.15;
        
    case{7} % Shocktube problem of G.A. Sod, JCP 27:1, 1978 
        fprintf('Shocktube problem of G.A. Sod, JCP 27:1, 1978 \n');
        p   = [1.0  0.1  ];
        u   = [0.75 0    ];
        rho = [1    0.125];
        tEnd = 0.10; cfl = 0.15;
        
    case{8} % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997 \n');
        p   = [3.528 0.571];
        u   = [0.698 0    ];
        rho = [0.445 0.5  ];
        tEnd = 0.08; cfl = 0.15; 
      
    case{9} % Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997 \n');
        p   = [10.333  1  ];
        u   = [ 0.92  3.55];
        rho = [ 3.857  1  ];
        tEnd = 0.09; cfl = 0.15;
        
    case{10} % Shocktube problem with supersonic zone
        fprintf('Shocktube problem with supersonic zone \n');
        p   = [1  0.02];
        u   = [0  0.00];
        rho = [1  0.02];
        tEnd = 0.162; cfl = 0.15; 

    case{11} % Contact discontinuity
        fprintf('Contact discontinuity \n');
        p   = [0.5 0.5];
        u   = [0.0 0.0];
        rho = [1.0 0.6];
        tEnd = 1; cfl = 0.15; 

    case{12} % Stationary shock
        fprintf('Stationary shock \n');
        p   = [ 1.0  0.1 ];
        u   = [-2.0 -2.0 ];
        rho = [ 1.0 0.125];
        tEnd = 0.1; cfl = 0.15;
        
	case{13} % Bogus case 1, Unknown
        fprintf('Case 7: Bogus case 1, Unknown \n');
        u  = [0         0        ];
        t  = [4.38385   8.972544 ];
        r  = [0.2253353 0.1204582];
        tEnd = 0.1; cfl = 0.15;
        
    case{14} % Manuel's Blastwave
        fprintf('Manuel`s Blastwave, IAM 2013 \n');
        p   = [4.00 0.10 1.00];
        u   = [1.00 0.00 -1.0];
        rho = [0.90 0.90 0.90];
        tEnd = 0.20; cfl = 0.15;
        
    otherwise 
        error('Available cases: 1~13');
        
end
% Compute Semiclassical ICs
if input ~= 13 
    E = p+(0.5).*rho.*u.^2; % Energy
    t = 4*E./rho-2*u.^2;    % Temperature
    r = rho./sqrt(pi*t);    % Fugacity
else
    % Do nothing
end

%% Load Selected case Initial condition:
% number of points required
 [k,nx] = size(x);

 switch input
     case{1,2,3,4,5,6,7,8,9,10,11,12,13} % L and R IC's
         % Preallocating
         r_0 = zeros(k,nx); u_0 = zeros(k,nx); t_0 = zeros(k,nx);
         p_0 = zeros(k,nx); rho_0 = zeros(k,nx); E_0 = zeros(k,nx);
         % Parameters of regions dimensions
         x_middle = (x(end)-x(1))/2;
         L = find(x<=x_middle);
         R = find(x>x_middle);
         
         % Initial Condition for our 2D domain
         % Fugacity
         r_0(L) = r(1); % region 1
         r_0(R) = r(2); % region 2
         % Velovity in x
         u_0(L) = u(1); % region 1
         u_0(R) = u(2); % region 2
         % temperature
         t_0(L) = t(1); % region 1
         t_0(R) = t(2); % region 2
         % pressure
         p_0(L) = p(1); % region 1
         p_0(R) = p(2); % region 2
         % density
         rho_0(L) = rho(1); % region 1
         rho_0(R) = rho(2); % region 2
         % Energy
         E_0(L) = E(1); % region 1
         E_0(R) = E(2); % region 2
         
     case{14} % 3-zones IC's
         % Preallocating
         r_0 = zeros(k,nx); u_0 = zeros(k,nx); t_0 = zeros(k,nx);
         p_0 = zeros(k,nx); rho_0 = zeros(k,nx); E_0 = zeros(k,nx);
         % Parameters of regions dimensions
         L = abs(x(end)-x(1));
         x_a = x(1); x_b = 0.3*L; x_c = 0.7*L; x_d = x(end);
         L = find(x>=x_a & x<=x_b);
         C = find(x> x_b & x<=x_c);
         R = find(x> x_c & x<=x_d);
         
         % Initial Condition for our 2D domain
         % Fugacity
         r_0(L) = r(1); % region 1
         r_0(C) = r(2); % region 2
         r_0(R) = r(3); % region 2
         % Velovity in x
         u_0(L) = u(1); % region 1
         u_0(C) = u(2); % region 2
         u_0(R) = u(3); % region 2
         % temperature
         t_0(L) = t(1); % region 1
         t_0(C) = t(2); % region 2
         t_0(R) = t(3); % region 2
         % pressure
         p_0(L) = p(1); % region 1
         p_0(C) = p(2); % region 2
         p_0(R) = p(3); % region 2
         % density
         rho_0(L) = rho(1); % region 1
         rho_0(C) = rho(2); % region 2
         rho_0(R) = rho(3); % region 2
         % Energy
         E_0(L) = E(1); % region 1
         E_0(C) = E(2); % region 2
         E_0(R) = E(3); % region 2
 end