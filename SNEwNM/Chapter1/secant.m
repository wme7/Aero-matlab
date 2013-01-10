function [x, hist] = secant(x, f, tola, tolr) 
%
%SECANT  This a secant code with no line search. 
%
% This code terminates on small relative-absolute errors
%
% function [X, HIST] = SECANT(X, F, TOLA, TOLR) 
%
% Inputs: X = initial iterate
%         F = function
%         TOLA = absolute error tolerance
%         TOLR = relative error tolerance
%
% Output: X = approximate result
%      HIST = array of iteration history, useful for tables and plots
%                The two columns are iteration number and residual.
%
%             If you leave the hist array  out by calling secant as
%             Z = SECANT(X, F, TOLA, TOLR) 
%             the storage for hist is not allocated and
%             the iteration history is not stored.
%
% Limit iterations to 100.
%
maxit=100;    
if nargout == 2
    histmp=zeros(maxit,2);
end
%
% Initialize
%
itc=0;
fc=feval(f,x);
%
% A secant code needs two iterations and the
% function values at each. I'll start it off with the
% old value xm = .99 x, if x ~= 0, or xm=.001, if x = 0,
% where x is the initial iterate.
%
if x~=0
    xm=x*.99;
else
    xm=.001;
end
fm=feval(f,xm);
%
tol=tolr*abs(fc)+tola;
%
%    Store iteration history?
%
if nargout == 2
    histmp(itc+1,:)=[itc, fc];
end
%
%    Main Secant loop.
%
while(abs(fc) > tol)
%
%    Iteration becoming unbounded?
%
    if(abs(x) > 1.d9)
        disp(' iterate too large');
        if nargout == 2
            histmp(itc+1,:)=[itc, fc];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
% Use a secant approximate derivative.
%
    del=x-xm;
    if del ~= 0
        fp=(fc-fm)/del;
        s=-fc/fp;
    else
        s=0;
    end
%
% Continue the secant iteration. Make sure you update xm and fm.
%
    fm=fc; xm=x;
    x=x+s;
    fc=feval(f,x);
    itc=itc+1;
%
% Too many iterations?
%
    if(itc > maxit)
        disp(' maxit reached');
        if nargout == 2
            histmp(itc+1,:)=[itc, fc];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
    if nargout == 2
        histmp(itc+1,:)=[itc,fc];
    end
end
if nargout == 2
    hist=histmp(1:itc+1,:);
end
