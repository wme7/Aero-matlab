function x = X_1d(z,n,j,E,theta,D) 
% Evaluate conservation of Energy
x = zeros(size(z));
switch theta
    case{1} % FD
        for i = 1:length(z) 
            x(i) = 2*E -D*(n/FD(z(i),D/2))^((2+D)/D) * FD(z(i),(D/2+1)) - j.^2/n;
        end
    case{-1} % BE
        for i = 1:length(z) 
            x(i) = 2*E -D*(n/BE(z(i),D/2))^((2+D)/D) * BE(z(i),(D/2+1)) - j.^2/n;
        end
    otherwise
        error('wrong statistic')
end