function F = stream3d(A,vector)
%% Stream function for any D2Qn Lattice
% There are only 4 basic movements and any other is a combination of the
% basic ones:

% Matlab equivalent solution:    
% F = circshift(A,[z,y,x])

%% Evaluate input data

v = length(vector);
    z = vector(1);
    y = vector(2);
    x = vector(3);

[l,m,n] = size(A);
    a = zeros(1,m,n);
    b = zeros(l,1,n);
    c = zeros(l,m,1);

%% Stream in XY plane
    
if x == 1       % x+ stream
    b = A(:,:,n);
    for i = n:-1:2;
        A(:,:,i) = A(:,:,i-1);
    end
    A(:,:,1) = b;
elseif x == -1  % x- stream
    b = A(:,:,1);
    for i = 1:n-1;
        A(:,:,i) = A(:,:,i+1);
    end
    A(:,:,n) = b;
else            % x=0 stream
    %do nothing
end

if y == 1       % y+ stream
    c = A(:,m,:);
    for j = m:-1:2;
        A(:,j,:) = A(:,j-1,:);
    end
    A(:,1,:) = c;
elseif y == -1  % y- stream
    c = A(:,1,:);
    for j = 1:m-1;
        A(:,j,:) = A(:,j+1,:);
    end
    A(:,m,:) = c;
else            % y=0 stream 
    %do nothing 
end

if z == 1       % z+ stream
    a = A(l,:,:);
    for k = l:-1:2;
        A(k,:,:) = A(k-1,:,:);
    end
    A(1,:,:) = a;
elseif z == -1  % z- stream
    a = A(1,:,:);
    for k = 1:l-1;
        A(k,:,:) = A(k+1,:,:);
    end
    A(l,:,:) = a;
else            % z=0 stream 
    %do nothing 
end

F = A;