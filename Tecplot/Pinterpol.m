function E = Pinterpol(x1,y1,x2,y2,ne)
%% Edge Point interpolation function.
delta_x = (x2-x1)/ne;
delta_y = (y2-y1)/ne;

    for i = 1:ne+1
        E(i,1) = x1 + (i-1)*delta_x;
        E(i,2) = y1 + (i-1)*delta_y;
    end
return