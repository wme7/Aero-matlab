function [a,b,c] = apply_DG_DOM(a,b,c,nv)
% Repeat arrays 'nv' times
    a = repmat(a,[1,1,nv]);     c = repmat(c,[1,1,nv]);
    b = repmat(b,[1,1,nv]);
    