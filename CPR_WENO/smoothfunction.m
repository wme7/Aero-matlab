% Smooth Measurement function
tic
term = zeros(1,size(ulp,2));
for j = tCells;
    
    for s = 1:K
        dpl = ulp(j-1,:);
        % Derivate 's' times
        for i = 1:s
            dpl = polyder(dpl);
        end
        integ = dx^(2*s-1)*polyint(conv(dpl,dpl));
        term(s) = polyval(integ,x(1+K,j-1))-polyval(integ,x(1,j-1));
    end
    Beta0 = sum(term);
    
    for s = 1:K
        dpl = ulp(j,:);
        % Derivate 's' times
        for i = 1:s
            dpl = polyder(dpl);
        end
        integ = dx^(2*s-1)*polyint(conv(dpl,dpl));
        term(s) = polyval(integ,x(1+K,j))-polyval(integ,x(1,j));
    end
    Beta1 = sum(term);
    
    for s = 1:K
        dpl = ulp(j+1,:);
        % Derivate 's' times
        for i = 1:s
            dpl = polyder(dpl);
        end
        integ = dx^(2*s-1)*polyint(conv(dpl,dpl));
        term(s) = polyval(integ,x(1+K,j+1))-polyval(integ,x(1,j+1));
    end
    Beta2 = sum(term);
end
toc
[Beta0;Beta1;Beta2]