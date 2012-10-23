function [r_0,u_0,t_0] = reimann_IC1d(x,input)
% Load the IC of a classical 2D Riemann Problem configuration.
% In the notation we take advantage of the matlab array notation as follows
%
% prop = [prop_left , prop_right]
%
% r = fugacity
% u = velocity in x direction
% t = temperature

%% Initial Physical Properties per case:
switch input
    case{1} % Configuration 1
        u = [0          0         ];
        t = [4.383850   8.972544  ];
        r = [0.2253353  0.1204582 ];
        
    case{2} % Configuration 2
        u = [0          0         ];
        t = [4.383850   8.972544  ];
        r = [0.2253353  0.1204582 ];
        
    case{3} % Configuration 3
        u = [0          0         ];
        t = [4.383850   8.972544  ];
        r = [0.2253353  0.1204582 ];
        
    case{4} % Configuration 4
        u = [0          0         ];
        t = [4.383850   8.972544  ];
        r = [0.2253353  0.1204582 ];
        
    otherwise 
        error('Available cases: 1, 2, 3 and 4');
        
end

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