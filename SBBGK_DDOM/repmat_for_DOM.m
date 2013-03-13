a = magic(4); % nx = ny = 4; nx*ny = 16
%a = reshape(a,1,3*3); % step not necessary!
a = reshape(a,1,1,4,4); % nv = 3 and
a = repmat(a,[3,3,1]);