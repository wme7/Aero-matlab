function [a,b,c] = apply_DOM(a,b,c,nv)
% Repeat arrays 'nv' times
    a = repmat(a,nv,1);     c = repmat(c,nv,1);
    b = repmat(b,nv,1);
    