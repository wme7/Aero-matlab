function [x, hist] = chordsol(x, f, tola, tolr)
%CHORDSOL  Chord method code for scalar equations.
%
% This code terminates on small relative-absolute errors
%
% function [X, HIST] = CHORDSOL(X, F, TOLA, TOLR)
%
% Inputs: X = initial iterate
%         F = function
%         TOLA = absolute error tolerance
%         TOLR = relative error tolerance
%
% Output: x=approximate result
%      HIST = array of iteration history, useful for tables and plots
%                The two columns are iteration number and residual.
%
%                If you leave this out by calling chordsol as
%                Z = CHORDSOL(X, F, TOLA, TOLR) 
%                the storage for hist is not allocated and
%                the iteration history is not stored.
%
% Uses: fordiff.m, Tim's numerical derivative
%
%
maxit=20;    % maximum iterations
alpha=1.d-4; 
small=1.d-8;
if nargout == 2
    histmp=zeros(maxit,2);
end
%
% Initialize
%

itc=0;
fc=feval(f,x);
fp=fordiff(f,x,fc);
tol=tolr*abs(fc)+tola;
%
%    Store iteration history?
%
if nargout == 2
    histmp(itc+1,:)=[itc, abs(fc)];
end
%
%    Main chord loop.
%
while(abs(fc) > tol)
%
%    Iteration becoming unbounded?
%
    if(abs(x) > 1.d9)
        disp(' iterate too large');
        if nargout == 2
            histmp(itc+1,:)=[itc,abs(fc)];
            hist=histmp(1:itc,:);
        end
    return
    end
%
%    Check for small relative derivatives.
%
    if(abs(fp) <= small*abs(fc))
        disp(' small derivative error');
    return;
    end
    s=-fc/fp;
    xt=x+s;
    ft=feval(f,xt);
    x=xt;
    itc=itc+1;
%
%    Too many iterations?
%
        if(itc > maxit)
            disp(' maxit reached');
            if nargout == 2
                histmp(itc+1,:)=[itc,abs(fc)];
                hist=histmp(1:itc,:);
            end
            return
        end
%
        fc=feval(f,x);
    if nargout == 2
        histmp(itc+1,:)=[itc,abs(fc)];
        end
    end
        if nargout == 2
        hist=histmp(1:itc+1,:);
        end
%
% A simple forward difference in MATLAB.
%
function fp = fordiff(f,x,fc)
        h=1.d-7;
%
%       Notice how f is evaluated using the feval command.
%
        fright=feval(f,x+h);
        fp=(fright-fc)/h;

