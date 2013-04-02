function dx = dX(z,n,theta,D) 
% Evaluate conservation of Energy
dx = zeros(size(z));
switch theta
    case{1} % FD
        for i = 1:length(z) 
            dx(i) = BE(z(i),1.5)/z(i)*n^(1.5)*BE(z(i),-0.5)*3/z(i)*BE(z(i),0.5)^(-4)-...
                1/z(i)*BE(z(i),0.5)^(-0.5)*n^3
        end
    case{-1} % BE
        for i = 1:length(z) 
            dx(i) = BE(z(i),1.5)/z(i)*n^(1.5)*BE(z(i),-0.5)*3/z(i)*BE(z(i),0.5)^(-4)-...
                1/z(i)*BE(z(i),0.5)^(-0.5)*n^3;
        end
    otherwise
        error('wrong statistic')
end