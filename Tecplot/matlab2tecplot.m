%% Write data from matlab to Tecplot 360
clear all; close all;

%% Generate data
x = 1:0.1:2; 
y = 0:0.2:1; 

[X,Y] = meshgrid(x,y);

Z = sin(y'*x);
[ny,nx] = size(Z);

a = 1; % zone or time step

%% Export to Tecplot
file = fopen('data1.plt','w');
% h1 gets the handel for the file "mesh1.tec".
% 'w' specifies that it will be written.
% similarly 'r' is for reading and 'a' for appending.
fprintf(file, 'TITLE = "Z = sin(YX)"\n');
fprintf(file, 'VARIABLES = "X", "Y", "Z"\n');


fprintf(file, 'ZONE T = "%1.4f"\n',a);
fprintf(file, 'I = %d, J = %d, K = 1, F = POINT\n\n', nx,ny);
for j = 1:ny
    for i = 1:nx
        %fprintf(file, '%f\t%f\t%f\n', x(i),y(j),Z(j,i));
        %or
        fprintf(file, '%f\t%f\t%f\n', X(1,i),Y(j,1),Z(j,i));
    end
end

fclose(file);