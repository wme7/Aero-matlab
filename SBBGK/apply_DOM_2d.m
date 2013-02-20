function [a,b,c,d] = apply_DOM_2d(a,b,c,d,nv,nx,ny)
% Matrix Array to Perform DOM operations
a = reshape(a,1,1,ny,nx); a = repmat(a,[nv,nv,1]);
b = reshape(b,1,1,ny,nx); b = repmat(b,[nv,nv,1]);
c = reshape(c,1,1,ny,nx); c = repmat(c,[nv,nv,1]);
d = reshape(d,1,1,ny,nx); d = repmat(d,[nv,nv,1]);