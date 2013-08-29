%% Just fot Testing
% symbolic manipulation
k = 3

% Solution points in the standard element [-1,1]
for i=1:k+1
    x0(i)= -cos((i-1)*pi/k);
end

x = sym('x');

% Lagrange polynomials
for i=1:k+1
    l(i)=x/x;
    for j=1:k+1
        if(i ~= j) 
            l(i)=l(i)*(x-x0(j))/(x0(i)-x0(j));
        end
    end
end
