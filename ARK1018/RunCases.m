%% Run Cases for ploting Data
clc; clear all; close all; 

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end

%% Number of cases
cases = 1:4;                    % Number of Cases to evaluate

%% Common Parameters
Nx          = 8;               % Number of elements
t_end       = 0.12;             % Final Time 
P_deg       = 4;                % Polinomial Degree
Pp          = P_deg+1;          % Polinomials Points
CFL         = 1/(2*P_deg+1);    % Courant Number 
RK_stages   = 4;                % Number of RK stages
plot_fig    = 0;                % {1}: plot while computing, {0}: no plot
theta       = 0;                % {-1} BE, {0} MB, {1} FD.

%% particular Paramerters
tau = [0.1 0.01 0.001 0.0001];  % Relaxation time per case

%% Initial Condition IC
ch=0;
while(ch==0)
fprintf ('Choose one of the following cases :- \n');
fprintf ('\n \t Case 1: Sods problem \n');
fprintf ('\t Case 2: Left running expansion and rightrunning "STRONG" shock \n');
fprintf ('\t Case 3: Left running shock and right runningexpansion \n');
fprintf ('\t Case 4: Double shock \n');
fprintf ('\t Case 5: Double expansion \n');
fprintf ('\t Case 6: Cavitation \n');
cas=input ('\nEnter a case no. <1-6>: ');
if cas==1
% Case 1:Left Expansion & right Shock
fprintf('Case 1:Sods problem \n');
rho1=1; rho4=0.125;
u1=0;   u4=0;
p1=1;   p4=0.1;
ch=1;
elseif cas==2
% Case 2:Strong Expansion & Shock
fprintf('Case 2:Strong ...Expansion & Shock \n');
rho1=3;  rho4=2;
u1=0;    u4=0;
p1=1000; p4=0.01;
ch=1;
elseif cas==3
% Case 3:Shock & Expansion
fprintf('Case 3:Shock & Expansion \n');
rho1=1; rho4=1;
u1=0;   u4=0;
p1=7;   p4=10;
ch=1;
elseif cas==4
% Case 4:Double Shock
fprintf('Case 4:Double Shock \n');
rho1=6; rho4=6;
u1=20;  u4=-6;
p1=450; p4=45;
ch=1;
elseif cas==5
% Case 5:Double Expansion
fprintf('Case 5:Double Expansion \n');
rho1=1; rho4=2.5;
u1=-2;  u4=2;
p1=40;  p4=40;
ch=1;
elseif cas==6
% Case 6:Cavitation
fprintf('Case 6:Cavitation \n');
rho1=1;  rho4=1;
u1=-20;  u4=20;
p1=0.40; p4=0.40;
ch=1;
else
fprintf ('Please enter an appropriate choice \n');
end % for case selection if-else loop
end % for case selection while loop
fprintf ('P1 = %f \n',p1);
fprintf ('P4 = %f \n',p4);
fprintf ('U1 = %f \n',u1);
fprintf ('U4 = %f \n',u4);
fprintf ('rho1 = %f \n',rho1);
fprintf ('rho4 = %f \n',rho4);

%% Compute Exact solution
[xt,u,rho,p,e] = Exact_Riemann(u1,p1,rho1,u4,p4,rho4);
    % plot exact solution
    if plot_fig == 1
        % Plotting instructions
        subplot (2,2,1)
        plot(xt,u);
        title('Plot of U v/s x/t');
        xlabel ('x/t');
        ylabel ('u');
        axis tight;
        subplot (2,2,2)
        plot(xt,rho);
        title('Plot of Density v/s x/t');
        xlabel ('x/t');
        ylabel ('rho');
        axis tight;
        subplot (2,2,3)
        plot(xt,p);
        title('Plot of Pressure v/s x/t');
        xlabel ('x/t');
        ylabel ('P');
        axis tight;
        subplot (2,2,4)
        plot(xt,e);
        title('Plot of Internal Energy v/s x/t');
        xlabel ('x/t');
        ylabel ('E');
        axis tight;
    end


%% Excecute SBBGK in Parallel using 4 processors
parfor i = cases
    [x(i,:),R(i,:),U(i,:),E(i,:),P(i,:),T(i,:),Z(i,:)] = ...
        DG_bgk1D_ERK(t_end,tau(i),Nx,P_deg,Pp,RK_stages,CFL,theta,plot_fig);
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