function B=colshift(A,k)
% Shifts the colums of A by k. So that the i-th column is shifted k(i) 
[m,n]=size(A);
% index vectors for rows and columns 
p = 1:m;
q = 1:n;
% index matrices for rows and columns
[P, Q] = ndgrid(p, q);
% create a matrix with the shift values
K = repmat(k(:), [1 m])';
% update the matrix with the row indexes
P =max(1,min(P-K,m));
% create matrix of linear indexes
ind = sub2ind([m n], P, Q);
% finally, create the output matrix
B = A(ind);
