%% Material Derivate Subroutine
% To Evaluate:  Du/Dt = du/dt + c du/dx
clear all; close all; clc;

dx = 0.01; x = 0:dx:1;
dt = 0.01; t = 0:dt:2;
c = 0.2;

%IC
u_0 = sin(2*pi*(x));

%load IC
u = u_0;

du = zeros(size(u));
for time = t
    subplot(1,2,1); 
    plot(x,u);
    
        
    u_old = u;
    u = sin(2*pi*(x + c*time));
        
    du(2:end-1) = u(3:end) - u(1:end-2);
    
    DuDt = (u - u_old)/dt + c*du/(2*dx);
    
    subplot(1,2,2); 
    plot(x,DuDt);
    
    
    drawnow
end
hold off