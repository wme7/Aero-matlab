%% Evaluating Fugacity Equations
% F_nu = Fermi function
% B_nu = Bose function
clear all; close all; clc;

% Fugacity range
z = linspace(0,1,40); 

% Spatial parameters
E = 0.5; j = 1.0; n = 1.0; 

for i = 1:length(z)
    
    % Evaluate Fermi Functions
    tic; X_1(i) = 2*E -(n/FD(z(i),0.5))^3*FD(z(i),1.5) - j.^2/n; t(1) = toc;
    tic; X_2(i) = 2*E - 2*(n/FD(z(i),1.0))^2*FD(z(i),2.0) - (j+j).^2/n; t(2) = toc;
    tic; X_3(i) = 2*E - 3*(n/FD(z(i),1.5))^(5/3)*FD(z(i),2.5) - (j+j+j).^2/n; t(3) = toc;
    
    % Evaluate Fermi Functions
    tic; Y_1(i) = 2*E -(n/BE(z(i),0.5))^3*BE(z(i),1.5) - j.^2./n; y(1) = toc;
    tic; y_2(i) = 2*E - 2*(n/BE(z(i),1.0))^2*BE(z(i),2.0) - (j+j).^2./n; y(2) = toc;
    tic; Y_3(i) = 2*E - 3*(n/BE(z(i),1.5))^(5/2)*BE(z(i),2.5) - (j+j+j).^2./n; y(3) = toc;
    
end

% Show in screen
disp(t');
disp(y');

%line = {'ks--','kd--','ko--','k*--'};
line = {'r*--','g*--','m*--','k*--'};

% Plot comparisons
range=[0,1,1e-2,0]; figure(1)
subplot(1,2,1); box on; 
plot(z,X_1,line{1}); hold on; plot(z,X_2,line{2}); plot(z,X_3,line{3});
hold off; legend('X_1{z}','X_2{z}','X_3{z}',4); title('for Fermi-Diract');
xlabel('Fugacity, z'); ylabel('X_a(z)'); 
subplot(1,2,2); box on; 
plot(z,Y_1,line{1}); hold on; plot(z,y_2,line{2}); plot(z,Y_3,line{3});
hold off; legend('X_1{z}','X_2{z}','X_3{z}',4); title('for Bose-Einstein'); 
xlabel('Fugacity, z'); ylabel('X_a(z)'); 

