function p = bisection(f,a,b,tol)

% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code and very handy!

if f(a)*f(b)>0
    disp('Wrong choice bro')
else
    p = (a + b)/2;
    err = abs(b-a);
    while err >= tol
        if f(a)*f(p)<0
            b = p;
        else
            a = p;
        end
        p = (a + b)/2;
        err = abs(b-a);
    end
end