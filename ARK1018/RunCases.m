%% Run Cases for ploting Data
clc; clear all; close all; 

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end

%% Initial Setup 
t_end = 0.12;
tau = [0.1 0.01 0.001 0.0001];
plot = 0; % {1}: plot while computing, {2}: no plot

%% Excecute Parallel
parfor i = 1:4
    [x(i,:),R(i,:),U(i,:)] = DG_bgk1D_ERK(t_end,tau(i),plot);
end
matlabpool close

%% Write results to tecplot
nx = length(x);
file = fopen('DG1D_results.plt','w');
% 'file' gets the handel for the file "case.plt".
% 'w' specifies that it will be written.
% similarly 'r' is for reading and 'a' for appending.

fprintf(file, 'TITLE = "DG_bgk1D_ERK results"\n');
fprintf(file, 'VARIABLES = "X" "R1" "R2" "R3" "R4" "U1" "U2" "U3" "U4"\n');
fprintf(file, 'ZONE T = "Final Time %0.2f"\n', t_end);
fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);

for i = 1:nx
    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', x(1,i),...
        R(1,i),R(2,i),R(3,i),R(4,i), ...
        U(1,i),U(2,i),U(3,i),U(4,i));
end

fclose(file);
%% if everything is ok then, tell me:
fprintf('File has been saved!\n')