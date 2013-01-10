function [x, hist] = newtsol(x, f, tola, tolr,jdiff) 
%NEWTSOL   Newton-Armijo code for scalar equations.
%
% This code terminates on small relative-absolute errors
%
%         [X, HIST] = NEWTSOL(X, F, TOLA, TOLR, JDIFF) 
%
% Inputs: X = initial iterate
%         F = function
%         TOLA = absolute error tolerance
%         TOLR = relative error tolerance
%         JDIFF = 0, analytic derivative 
%                    syntax: [f,fp]=f(x)
%         JDIFF = 1, forward difference derivative 
%
% Output: x=approximate result
%      HIST = array of iteration history, useful for tables and plots
%            The four columns are iteration number, residual, 
%            number of step size reductions done in the line
%            search, and the Newton iterate itself
%
%
%             If you leave out the hist array by calling newtsol as
%             Z = NEWTSOL(X, F, TOLA, TOLR) 
%                then storage for the history array is not allocated and
%                the iteration history is not stored.
%
if nargin == 4
    jdiff = 1;
end
maxit=100;    % maximum iterations
maxitls=20;   % maximum iterations inside the line search
alpha=1.d-4; 
small=1.d-8;
%
% Initialize
%
itc=0;
fc=feval(f,x);
tol=tolr*abs(fc)+tola;
%
% Store iteration history?
%
if nargout == 2
    histmp=zeros(maxit,4);
    histmp(itc+1,:)=[itc, abs(fc), 0, x];
end
%
%    Main Newton loop
%
while(abs(fc) > tol)
%
%    Iteration becoming unbounded?
%
    if(abs(x) > 1.d9)
        disp(' iterate too large');
        if nargout == 2
            histmp(itc+1,:)=[itc, abs(fc), 0, x];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
    lambda=1;
    iarm=0;           % line search iteration counter
    if jdiff==1
        fp=fordiff(f,x,fc);% use a numerical derivative
    else
        [fv,fp]=feval(f,x);
    end
%
%    Check for small relative derivatives.
%
    if(abs(fp) <= small*abs(fc))
        disp(' small derivative error');
        return;
    end
    s=-fc/fp;
    xt=x+lambda*s;
    ft=feval(f,xt);
%
%    Main line search loop. 
%
    while(abs(ft) >= (1 - alpha*lambda)*abs(fc) + 1.d-12)
        lambda=lambda/2;
        xt=x+lambda*s;
        ft=feval(f,xt);
        iarm=iarm+1;
%
%    Are you spending too much time in the line search?    
%
        if(iarm > maxitls)
            disp(' line search failure');
            if nargout == 2
                histmp(itc+1,:)=[itc, abs(fc), 0, x];
                hist=histmp(1:itc+1,:);
            end
            return
        end
    end
%
%    Step accepted, continue the Newton iteration.
%
    x=xt;
    itc=itc+1;
%
%    Too many iterations?
%
    if(itc > maxit)
        disp(' maxit reached');
        if nargout == 2
            histmp(itc+1,:)=[itc, abs(fc), 0, x];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
    fc=feval(f,x);
    if nargout == 2
       histmp(itc+1,:)=[itc,abs(fc),iarm,x];
    end
end
%
%   Fix up the history array.
%
if nargout == 2
    hist=histmp(1:itc+1,:);
end
%
% A simple forward difference in MATLAB.
%
function fp = fordiff(f,x,fc)
h=1.d-7;
%
%       Notice how f is evaluated using the feval command
%
fright=feval(f,x+h);
fp=(fright-fc)/h;

