%% FLow over a Rectangle
% This algorithm is based on the journal paper: 
% An accuratee Direct Solver for Semiclassical Botlzmann-BGK Equation using 
% high order Weno Scheme and Discrete Ordinate Method. by Bagus Putra
% Muljadi

% Manuel A. Diaz
% National Taiwan University
% Institute of Applied Mechanics
% Wednesday 14th, 2011.
clear;
clc;
close all;

%% Build Grid and Length of the Domain
% We initialy wish to use an 100x100 2D Grid
n=101;              %Number of grid points
L=10;               %Length of domain
h=L/(n-1);          %Spatial step size
x=0:h:L;            %X axis values

%% Base Constants
CFL=0.34;           %CFL number for stability
t_final=3.9e-3;     %Final time
gamma=1.4;          %Ratio of specific heats for ideal di-atomic gas

%% Select Solution Method
i_solver=0;     %Euler=0, N-S=1
i_exim=1;       %Explicit=0,    Implicit=1
i_method=0;     %Upwind & TVD=0. weno=1

%% Initial Conditions
p_l=1e5;            %Pressure in left side of shock tube at t=0
p_r=1e3;            %Pressure in right side of shock tuve at t=0
rho_l=1;            %Density at left side of shock tube at t=0
rho_r=0.01;         %Density at right side of shock tube at t=0
u_l=0;              %Velocity in left side of shock tube at t=0
u_r=0;              %Velocity in right side of shock tube at t=0

