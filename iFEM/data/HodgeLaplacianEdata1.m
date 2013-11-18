function pde = HodgeLaplacianEdata1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.f = @f;
pde.exactu = @exactu;
pde.exactsigma = @exactsigma;
pde.gu = @exactu;
pde.gsigma = @exactsigma;
pde.gun = @gun;
pde.grotu = 0;

    function s = f(p)
    s = exactu(p);
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2);
    s = [cos(x), -sin(y)];
    end

    function s = exactsigma(p)
    x = p(:,1); y = p(:,2);
    s = sin(x)+cos(y);   
    end

    function f = gun(p) % for unit square [0,1]^2
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    u = exactu(p);
    leftbd = (abs(x)<eps);  % n = (-1,0); 
    f(leftbd) = - u(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0); 
    f(rightbd) = u(rightbd,1);
    topbd = (abs(y-1)<eps);   % n = (0,1)
    f(topbd) = u(topbd,2);
    bottombd = (abs(y)<eps);% n = (0,-1)
    f(bottombd) = - u(bottombd,2);
    end

end