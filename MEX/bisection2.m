function p = bisection2(a,b,tol) 


% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code and very handy!

if myfunc(a)*myfunc(b)>0
    error('Wrong choice bro');
else
    p = (a + b)/2;
    err = abs(b-a);
    while err >= tol
        if myfunc(a)*myfunc(p)<0
            b = p;
        else
            a = p;
        end
        p = (a + b)/2;
        err = abs(b-a);
    end
end