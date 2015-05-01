% collision test

% ball 1
x1 = 0;
y1 = 0;
u1 = 0;
v1 = 0;
d1 = 5;
m1 = 1;

% ball 2
x2 = 20;
y2 = 20;
%u2 = 
%v2 = 
d2 = 5;
m2 = 1;

dist12 = sqrt( (x1-x2)^2 + (y1-y2)^2 );

% compute momentum before collision
px = m1*u1 + m2*u2;
py = m1*v1 + m2*v2;

if dist12 == (d1+d2)/2
    disp('collision')
    % compute momentum after collision
    
else
    disp('far away')
    % do nothing
end

h=scatter([x1,x2],[y1,y2]); 
axis([-0.5,20.5,-0.5,20.5]);