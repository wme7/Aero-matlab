%% Vectorization
% 
% We summarize common tricks of vectorization used in iFEM. We assume
% readers know basic tricks on the performance improvment; search
% *Techniques for Improving Performance* in Matlab help window. Scroll
% down and read *Simple Example of Vectorizing* and *Advanced Example of
% Vectorizing*

%% Linear indexing
%
% For a vector, to access a range of elements, you can use a logical array.
% A typical example:

a = rand(10,1);
idx = (a>0.5);
apositive = a(idx);
display(apositive);
%%
% For a matrix, unfortunately if we use two indices arrays in the subscript,
% it will extract a submatrix.

A = magic(6);
display(A);
i = [1 2 3];
j = [3 4 5];
display(A(i,j));
%%
% What we want is indeed A(1,3), A(2,4) and A(3,5) not the submatrix
% A(1:3,3:5).
%
% The solution is to transfer subscript to linear indexing. Every matrix
% can be accessed with a single subscript, |A(k)|, since the matrix is
% stored as a long array stacked by columns. Search *Linear Indexing* in
% Matlab help window.
%
% The function |sub2ind| is designed for this purpose. 
idx = sub2ind(size(A),i,j);
display(A(idx));
%%
% The sub2ind only works for dense matrices not sparse matrices since the
% indexing of sparse matrices is much more involved. For sparse matrices,
% use spsub2ind to extract the values.
n = 5;
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
display(full(A));
i = [1 2 4];
j = [3 1 4];
Aij = spsub2ind(A,i,j);
display(i); display(j);
display(Aij);
%%
% Again A(i,j) will extract a sub-matrix. 
display(full(A(i,j)));