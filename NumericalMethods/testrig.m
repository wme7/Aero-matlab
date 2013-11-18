
% list of number of nodes to use for each run
Ms = [20 40 80 160 320];

for cnt = 1:5
    M = Ms(cnt);
    
    centraldifference2
    
    saverrors{cnt} = finalerror;
    savex{cnt} = x;
end

% plot errors
semilogy(savex{1}, saverrors{1}, 'r-*');
hold on;
semilogy(savex{2}, saverrors{2}, 'k-s');
semilogy(savex{3}, saverrors{3}, 'g-o');
semilogy(savex{4}, saverrors{4}, 'b-d');
semilogy(savex{5}, saverrors{5}, 'm-+');
hold off;

legend('M=20', 'M=40', 'M=80', 'M=160', 'M=320');
xlabel('x', 'FontSize', 16);
ylabel('|u_{exact}-u_{numerical}|', 'FontSize', 16);

% compute maximum absolute error

for cnt=1:5
    error(cnt) = max(abs(saverrors{cnt}));
end
dx = 2./Ms;
figure; loglog(dx, error, 'r-*'); 
hold on; loglog(dx, dx.^2, 'k-'); hold off;
xlabel('dx', 'FontSize', 16);
ylabel('max(|u_{exact}-u_{numerical}|)', 'FontSize', 16);
title('Maximum error at t=8');
legend('error', 'reference dx^2');
