%% Mesh Generator for Tecplot
% Generating a 2D mesh for Tecplot for an axisymmetric diffuser.
% by Manuel Diaz, 2012.09.15.

clear all; close all; 
%% Parameters
L0 = 0; L1 = 10 + L0; L2 = 6 + L1 + L0;
H0 = 0; H1 = 0.5;
theta = deg2rad(5); % theta = 5 degress converted to radians
H2 = H1 + L2*tan(theta);

%% Number of elements of edges (ne)
ne_x1 = 25; ne_x2 = 40; ne_y = 40;

%% End Points for the rectangular zones:
Z1=[ L0,H1 ; L1,H1 ; L1,H0 ; L0,H0 ]; % Zone1
Z2=[ L1,H1 ; L2,H2 ; L2,H0 ; L1,H0 ]; % Zone2

%% Zone 1: 
% There are 4 edge lines, each having several edge points.
% We use the Function "Pinterpol" to create the edge points for zone 1.

E1 = Pinterpol(Z1(1,1),Z1(1,2),Z1(2,1),Z1(2,2),ne_x1);
E2 = Pinterpol(Z1(4,1),Z1(4,2),Z1(3,1),Z1(3,2),ne_x1);

% Now generating the mesh for Zone 1
M1 = zeros(ne_x1+1,ne_y+1,2);
for i = 1:(ne_x1+1);
    V = Pinterpol(E1(i,1),E1(i,2),E2(i,1),E2(i,2),ne_y);
    M1(i,1:(ne_y+1),1:2) = V; % mesh1
end

%% Zone 2:
% Similarly for zone 2
E1 = Pinterpol(Z2(1,1),Z2(1,2),Z2(2,1),Z2(2,2),ne_x2);
E2 = Pinterpol(Z2(4,1),Z2(4,2),Z2(3,1),Z2(3,2),ne_x2);

% generating the mesh for Zone 2
M2 = zeros(ne_x2+1,ne_y+1,2);
for i = 1:(ne_x2+1);
    V = Pinterpol(E1(i,1),E1(i,2),E2(i,1),E2(i,2),ne_y);
    M2(i,1:(ne_y+1),1:2) = V; % mesh2
end

%% To Tecplot
% Writing Output data for Tecplot format:

file = fopen('mesh1.tec','w');
% h1 gets the handel for the file "mesh1.tec".
% 'w' specifies that it will be written.
% similarly 'r' is for reading and 'a' for appending.

fprintf(file, 'TITLE = "Mesh from Matlab"\n');
fprintf(file, 'VARIABLES = "X" "Y"\n');
fprintf(file, 'ZONE T = "Inlet zone"\n');
fprintf(file, 'I = %d, J = %d, K = 1, F = POINT\n\n', ne_x1+1,ne_y+1);

for j = 1:(ne_y+1)
    for i = 1:(ne_x1+1)
        fprintf(file, '%f\t%f\n', M1(i,j,1),M1(i,j,2));
    end
end

fprintf(file, 'ZONE T="Diffuser zone"\n');
fprintf(file, 'I = %d, J=%d, k=1, F=POINT\n\n', ne_x2+1,ne_y+1);

for j = 1:(ne_y+1)
    for i = 1:(ne_x2+1)
        fprintf(file, '%f\t%f\n', M2(i,j,1),M2(i,j,2));
    end
end

fclose(file);