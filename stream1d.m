function F = stream1d(A,vector)
%% Stream function for any D1Qn Lattice
% There are only 2 basic movements and any other is a combination of the
% basic ones:

% Matlab equivalent solution:    
% F = circshift(A,[y,x])

%% Evaluate input data

v = length(vector);

    x = vector(2);
    y = vector(1);

n = length(A);

    b = zeros(1,n);

%% Stream in XY plane
    
if x == 1       % x+ stream
    b = A(:,n);
    for i = n:-1:2;
        A(:,i) = A(:,i-1);
    end
    A(:,1) = b;
elseif x == -1  % x- stream
    b = A(:,1);
    for i = 1:n-1;
        A(:,i) = A(:,i+1);
    end
    A(:,n) = b;
else            % x=0 stream
    %do nothing
end

F = A;