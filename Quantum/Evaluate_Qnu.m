%% Evaluating Quantum functions, Q_nu
% F_nu = Fermi function
% B_nu = Bose function
clear all; close all; clc;

% Fugacity range
z = linspace(0,1,40); 

for i = 1:length(z)
    
    % Evaluate Fermi Functions
    tic; FD_12(i) = FD(z(i),0.5); t(1) = toc;
    tic; FD_32(i) = FD(z(i),1.5); t(2) = toc;
    tic; FD_52(i) = FD(z(i),2.5); t(3) = toc;
    
    % Evaluate Fermi Functions
    tic; BE_12(i) = BE(z(i),0.5); y(1) = toc;
    tic; BE_32(i) = BE(z(i),1.5); y(2) = toc;
    tic; BE_52(i) = BE(z(i),2.5); y(3) = toc;
    
end

% Show in screen
disp(t');
disp(y');

line = {'ks--','kd--','ko--','k*--'};
%line = {'r*--','g*--','m*--','k*--'};

% Plot comparisons
figure(1)
subplot(2,2,1); box on; xlabel('Fugacity, z'); ylabel('FD(z)')
hold on; plot(z,FD_12,line{1}); plot(z,FD_32,line{2}); plot(z,FD_52,line{3});
hold off; legend('FD_{1/2}','FD_{3/2}','FD_{5/2}',2);
subplot(2,2,2); box on; xlabel('Fugacity, z'); ylabel('BE(z)')
hold on; plot(z,BE_12,line{1}); plot(z,BE_32,line{2}); plot(z,BE_52,line{3});
hold off; legend('BE_{1/2}','BE_{3/2}','BE_{5/2}',2);
subplot(2,2,3); box on; xlabel('Fugacity, z'); ylabel('FD_a(z)/FD_b(z)')
hold on; plot(z,FD_32./FD_12,line{1}); plot(z,FD_52./FD_32,line{2}); plot(z,FD_52./FD_12,line{3});
hold off; legend('FD_{3/2}/FD_{1/2}','FD_{5/2}/FD_{3/2}','FD_{5/2}/FD_{1/2}',2);
subplot(2,2,4); box on; xlabel('Fugacity, z'); ylabel('BE_a(z)/BE_b(z)')
hold on; plot(z,BE_32./BE_12,line{1}); plot(z,BE_52./BE_32,line{2}); plot(z,BE_52./BE_12,line{3});
hold off; legend('BE_{3/2}/BE_{1/2}','BE_{5/2}/BE_{3/2}','BE_{5/2}/BE_{1/2}',3);
