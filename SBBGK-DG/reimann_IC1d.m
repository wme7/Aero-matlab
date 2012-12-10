function [r_0,u_0,t_0] = reimann_IC1d(x,input)
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
% E = Enerty
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
        p   = [1    0.1  ];
        u   = [0.75 0    ];
        rho = [1    0.125];
        
    case{2} % Configuration 2, Left Expansion and right strong shock
        p   = [1000 0.1  ];
        u   = [0    0    ];
        rho = [3    2    ];
        
    case{3} % Configuration 3, Right Expansion and left strong shock
        p   = [7    10   ];
        u   = [0    0    ];
        rho = [1    1    ];
        
    case{4} % Configuration 4, Double Shock
        p   = [450  45   ];
        u   = [20   -6   ];
        rho = [6    6    ];
        
    case{5} % Configuration 5, Double Expansion
        p   = [40   40   ];
        u   = [-2   2    ];
        rho = [1    2.5  ];

    case{6} % Configuration 6, Cavitation
        p   = [0.4  0.4  ];
        u   = [-20  20   ];
        rho = [1    1    ];
        
    case{7} % Bogus case 1
        % Case 1
        % UL  = 0.; UR  = 0.;
        % TL  = 4.38385; TR  = 8.972544;
        % ZL  = 0.2253353; ZR  = 0.1204582;
        
    case{8} % Bogus case 2
        % Case 2
        % UL  = 0.; UR  = 0.;
        % TL  = 4.38385; TR  = 8.972544;
        % ZL  = 0.2253353; ZR  = 0.1204582;

    otherwise 
        error('Available cases: 1, 2, 3, 4, 5 and 6');
        
end
% Compute Semiclassical ICs

% ET=PL+0.5*RL*UL^2;
% TL=4*ET/RL-2*UL^2;
% ZL=RL/sqrt(pi*TL);

E = p+(0.5).*rho.*u.^2; % Energy
t = 4*E./rho-2*u.^2;    % Temperature
r = rho./sqrt(pi*t);    % Fugacity

%% Load Selected case Initial condition:
% number of points required
 nx = length(x);

% Preallocating
r_0 = zeros(1,nx); 
u_0 = zeros(1,nx); 
t_0 = zeros(1,nx); 

% Parameters of regions dimensions
x_middle = ceil(nx/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;

% Initial Condition for our 2D domain
% Fugacity
r_0(l_1) = r(1); % region 1
r_0(l_2) = r(2); % region 2
% Velovity in x
u_0(l_1) = u(1); % region 1
u_0(l_2) = u(2); % region 2
% temperature
t_0(l_1) = t(1); % region 1
t_0(l_2) = t(2); % region 2