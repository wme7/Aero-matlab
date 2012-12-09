%% Run Cases for ploting Data
clc; clear all; close all; 

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end

%% Number of cases
cases = 1:4;                    % Number of Cases to evaluate

%% Common Parameters
Nx          = 32;               % Number of elements
t_end       = 0.12;             % Final Time 
P_deg       = 4;                % Polinomial Degree
Pp          = P_deg+1;          % Polinomials Points
CFL         = 1/(2*P_deg+1);    % Courant Number 
RK_stages   = 4;                % Number of RK stages
plot_fig    = 0;                % {1}: plot while computing, {0}: no plot
theta       = 0;                % {-1} BE, {0} MB, {1} FD.
iV          = 80;               % Space Velocity Points

%% particular Paramerters
tau = [0.1 0.01 0.001 0.0001];  % Relaxation time per case

%% Excecute SBBGK in Parallel using 4 processors

parfor i = cases
    [x(i,:),R(i,:),U(i,:),E(i,:),P(i,:),T(i,:),Z(i,:)] = ...
        DG_bgk1D_ERK_func(t_end,tau(i),Nx,P_deg,Pp,RK_stages,CFL,theta,plot_fig,iV);
end
matlabpool close

%% Create Distinctive ID names for each result case:
IDn1 = 0; IDn2 = 0; IDn3 = 0; IDn4 = 0;
ID1 = 0;  ID2 = 0;  ID3 = 0;  ID4 = 0;
for j = cases
    % Simulation Name
    name      = 'DG1D';
    % Statistic Used
    switch theta
        case {-1}
            statistic = 'BE';
        case {0}
            statistic = 'MB';
        case {1}
            statistic = 'FD';
        otherwise
            error('not a valid theta value')
    end
    % Elements used
    elements  = num2str(Nx);
    % Polinomial Degree
    P_degree = num2str(P_deg);
    % RK stages
    RKs = num2str(RK_stages);
    % Relaxation time
    w = 1/tau(j);
    w = num2str(w);
    % Number of velocity points
    
    % Tecplot format
    f = '.plt';
    % ID
    IDn = [name,statistic,'X',elements,' ','P',P_degree,'RK',RKs,'w',w,f];
    IDn4 = IDn3; IDn3 = IDn2; IDn2 = IDn1; IDn1 = IDn;
    ID = [name,statistic,'X',elements,'-','P',P_degree,'RK',RKs,'w',w];
    ID4 = ID3; ID3 = ID2; ID2 = ID1; ID1 = ID; 
end
IDn = {IDn4,IDn3,IDn2,IDn1};
ID = {ID4,ID3,ID2,ID1};

%% Write results to tecplot

for i = cases
    nx = length(x);
    % Open file
    file = fopen(IDn{i},'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
    
    fprintf(file, 'TITLE = "%s"\n',ID{i});
    fprintf(file, 'VARIABLES = "x" "Density" "Velocity in x" "Energy" "Pressure" "Temperature" "Fugacity"\n');
    fprintf(file, 'ZONE T = "Final Time %0.2f"\n', t_end);
    fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);
    
    for j = 1:nx
        fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
            x(i,j),R(i,j),U(i,j),E(i,j),P(i,j),T(i,j),Z(i,j));
    end
    % Close file
    fclose(file);
end
%% if everything is ok then, tell me:
fprintf('All Results have been succesfully saved!\n')