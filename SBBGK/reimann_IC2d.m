function [r_0,u_0,v_0,t_0] = reimann_IC2d(x,y,input)
% Load the IC of a classical 2D Riemann Problem configuration. 
% In the notation we take advantage of the matlab array notation as follows
%
% prop = [prop_reg1 , prop_reg2 , prop_reg3 , prop_reg4]
%
% r = fugacity
% u = velocity in x direction
% v = velocity in y direction
% t = temperature

%% Initial Physical Properties per case:
switch input
    case{5} % Configuration 5
        r = [ 0.142  0.4253 0.142  0.6635 ];
        u = [-0.75  -0.75   0.75   0.75   ];
        v = [-0.5    0.5    0.5   -0.5    ];
        t = [ 2.078  1.1494 2.078  0.87685];
        
    case{13} % Configuration 13
        r = [0.142   0.4253 0.4448 0.151  ];
        u = [0       0      0      0      ];
        v = [-0.3    0.3    0.697  0.26254];
        t = [2.0782  1.1494 0.7083 1.273  ];
        
    case{17} % Configuration 13
        r = [ 0.142011  0.142011 0.38035  0.142  ];
        u = [ 0.1       0.0      0.0      0.0    ];
        v = [-0.4      -0.3      0.1264  -0.97905];
        t = [ 2.07823   1.24705  0.77464  1.31443];
    otherwise 
        error('Available cases: 5, 13 and 17');
        
end

%% Load Selected case Initial condition:
% Parameters
[ny nx] = size(x);
[ny nx] = size(y);


% Preallocate u_0,
r_0 = zeros(ny,nx);
u_0 = zeros(ny,nx);
v_0 = zeros(ny,nx);
t_0 = zeros(ny,nx);

% Parameters of regions dimensions
x_middle = ceil(nx/2);
y_middle = ceil(ny/2);
l_1 = 1:x_middle; l_2 = x_middle+1:nx;
h_1 = 1:y_middle; h_2 = y_middle+1:ny;

% Initial Condition for our 2D domain
% Fugacity
r_0(h_1,l_1) = r(3); % region 1
r_0(h_1,l_2) = r(4); % region 2
r_0(h_2,l_1) = r(2); % region 3
r_0(h_2,l_2) = r(1); % region 4
% velocity in x
u_0(h_1,l_1) = u(3); % region 1
u_0(h_1,l_2) = u(4); % region 2
u_0(h_2,l_1) = u(2); % region 3
u_0(h_2,l_2) = u(1); % region 4
% velocity in y
v_0(h_1,l_1) = v(3); % region 1
v_0(h_1,l_2) = v(4); % region 2
v_0(h_2,l_1) = v(2); % region 3
v_0(h_2,l_2) = v(1); % region 4
% temperature
t_0(h_1,l_1) = t(3); % region 1
t_0(h_1,l_2) = t(4); % region 2
t_0(h_2,l_1) = t(2); % region 3
t_0(h_2,l_2) = t(1); % region 4
