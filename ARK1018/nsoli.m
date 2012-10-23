function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol, parms)
% NSOLI  Newton-Krylov solver, globally convergent 
%        solver for f(x) = 0
%
% Inexact-Newton-Armijo iteration
%
% Eisenstat-Walker forcing term
%
% Parabolic line search via three point interpolation.
%
% C. T. Kelley, April 27, 2001
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%        function = f
%        tol = [atol, rtol] relative/absolute
%            error tolerances for the nonlinear iteration
%        parms = [maxit, maxitl, etamax, lmeth, restart_limit]
%            maxit = maxmium number of nonlinear iterations
%                default = 40
%            maxitl = maximum number of inner iterations before restart
%                in GMRES(m), m = maxitl 
%                default = 40
%                
%                For iterative methods other than GMRES(m) maxitl
%                is the upper bound on linear iterations.
%
%            |etamax| = Maximum error tolerance for residual in inner
%                iteration. The inner iteration terminates
%                when the relative linear residual is
%                smaller than eta*| F(x_c) |. eta is determined
%                by the modified Eisenstat-Walker formula if etamax > 0.
%                If etamax < 0, then eta = |etamax| for the entire
%                iteration.
%                default: etamax = .9
%
%            lmeth = choice of linear iterative method
%                    1 (GMRES), 2 GMRES(m), 
%                    3 (BICGSTAB), 4 (TFQMR)
%                 default = 1 (GMRES, no restarts)
%
%            restart_limit = max number of restarts for GMRES if
%                    lmeth = 2
%                  default = 20
%
% output:
%        sol = solution
%        it_hist(maxit,3) = l2 norms of nonlinear residuals
%            for the iteration, number of function evaluations,
%            and number of steplength reductions
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
%
% internal parameters:
%       debug = turns on/off iteration statistics display as
%               the iteration progresses
%
%       alpha = 1.d-4, parameter to measure sufficient decrease
%
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%
%       maxarm = 20, maximum number of steplength reductions before
%                    failure is reported
%
%
%
%
% Set the debug parameter; 1 turns display on, otherwise off.
%
debug = 0;
%
% Set internal parameters.
%
alpha = 1.d-4; sigma0 = .1; sigma1 = .5; maxarm = 20; gamma = .9;
%
% Initialize it_hist, ierr, x_hist, and set the default values of
% those iteration parameters which are optional inputs.
%
ierr = 0; maxit = 40; lmaxit = 40; etamax = .9; it_histx = zeros(maxit,3);
lmeth = 1; restart_limit = 20;
if nargout == 4, x_hist = x; end
%
% Initialize parameters for the iterative methods.
% Check for optional inputs.
%
gmparms = [abs(etamax), lmaxit];
if nargin == 4
    maxit = parms(1); lmaxit = parms(2); etamax = parms(3);
    it_histx = zeros(maxit,3);
    gmparms = [abs(etamax), lmaxit];
    if length(parms)>= 4
       lmeth = parms(4);
    end
    if length(parms) == 5
       gmparms = [abs(etamax), lmaxit, parms(5), 1];
    end
end
%
rtol = tol(2); atol = tol(1); n = length(x); fnrm = 1; itc = 0;
%
% Evaluate f at the initial iterate,and
% compute the stop tolerance.
%
f0 = feval(f,x);
fnrm = norm(f0);
it_histx(itc+1,1) = fnrm; it_histx(itc+1,2) = 0; it_histx(itc+1,3) = 0;
fnrmo = 1;
stop_tol = atol + rtol*fnrm;
outstat(itc+1, :) = [itc fnrm 0 0 0];
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
%
% Keep track of the ratio (rat = fnrm/frnmo)
% of successive residual norms and 
% the iteration counter (itc).
%
    rat = fnrm/fnrmo;
    fnrmo = fnrm; 
    itc = itc+1;
    [step, errstep, inner_it_count,inner_f_evals] = ...
         dkrylov(f0, f, x, gmparms, lmeth);

%
%   The line search starts here.
%
    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    xt = x + lambda*step;
    ft = feval(f,xt);
    nft = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
    while nft >= (1 - alpha*lambda) * nf0;
%
%   Apply the three point parabolic model.
%
        if iarm == 0
            lambda = sigma1*lambda; 
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm); 
        end
%
% Update x; keep the books on lambda.
%
        xt = x+lambda*step;
        lamm = lamc;
        lamc = lambda;
%
% Keep the books on the function norms.
%
        ft = feval(f,xt);
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;
        iarm = iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            ierr = 2;
            disp(outstat)
            it_hist = it_histx(1:itc+1,:);
        if nargout == 4, x_hist = [x_hist,x]; end
            sol = xold;
            return;
        end
    end
    x = xt;
    f0 = ft;
%
%   End of line search.
%
    if nargout == 4, x_hist = [x_hist,x]; end
    fnrm = norm(f0);
    it_histx(itc+1,1) = fnrm; 
%
%   How many function evaluations did this iteration require?
%
    it_histx(itc+1,2) = it_histx(itc,2)+inner_f_evals+iarm+1;
    if itc == 1, it_histx(itc+1,2) = it_histx(itc+1,2)+1; end;
    it_histx(itc+1,3) = iarm;
%
    rat = fnrm/fnrmo;
%
%   Adjust eta as per Eisenstat-Walker.
%
    if etamax > 0
        etaold = gmparms(1);
        etanew = gamma*rat*rat;
        if gamma*etaold*etaold > .1
            etanew = max(etanew,gamma*etaold*etaold);
        end
        gmparms(1) = min([etanew,etamax]);
        gmparms(1) = max(gmparms(1),.5*stop_tol/fnrm);
    end
%
    outstat(itc+1, :) = [itc fnrm inner_it_count rat iarm];
%
end
sol = x;
it_hist = it_histx(1:itc+1,:);
if debug == 1
    disp(outstat)
    it_hist = it_histx(1:itc+1,:);
end
%
% on failure, set the error flag
%
if fnrm > stop_tol, ierr = 1; end
%
%
function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

%
% set internal parameters
%
sigma0 = .1; sigma1 = .5;
%
% compute coefficients of interpolation polynomial
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda
%
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end
%
%
%
function [step, errstep, total_iters, f_evals] = ...
    dkrylov(f0, f, x, params, lmeth)
% Krylov linear equation solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [step, errstep, total_iters, f_evals] 
%                              = dkrylov(f0, f, x, params, lmeth)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         x = current point
%         params = vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%              params(3) = max number of restarts for GMRES(m)
%              params(4) (Optional) = reorthogonalization method in GMRES
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         lmeth = method choice
%              1 GMRES without restarts (default)
%              2 GMRES(m), m = params(2) and the maximum number
%                   of restarts is params(3) 
%              3 Bi-CGSTAB
%              4 TFQMR
%
% Output: x = solution
%         errstep = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
%

%
% initialization
%
lmaxit = params(2);
restart_limit = 20;
if length(params) >= 3
    restart_limit = params(3);
end
if lmeth == 1, restart_limit = 0; end
if length(params) == 3
%
% default reorthogonalization
%
     gmparms = [params(1), params(2), 1];
elseif length(params) == 4
%
% reorthogonalization method is params(4)
%
     gmparms = [params(1), params(2), params(4)];
else
     gmparms = [params(1), params(2)];
end
%
% linear iterative methods
%
if lmeth == 1 | lmeth == 2  % GMRES or GMRES(m) 
%
% compute the step using a GMRES routine especially designed
% for this purpose
%
    [step, errstep, total_iters] = dgmres(f0, f, x, gmparms);
    kinn = 0;
%
%   restart at most restart_limit times
%
    while total_iters == lmaxit & ...
          errstep(total_iters) > gmparms(1)*norm(f0) & ...
          kinn < restart_limit
        kinn = kinn+1;
        [step, errstep, total_iters] = dgmres(f0, f, x, gmparms,step);
    end
    total_iters = total_iters+kinn*lmaxit;
    f_evals = total_iters+kinn;
%
% Bi-CGSTAB
%
elseif lmeth == 3
    [step, errstep, total_iters] = dcgstab(f0, f, x, gmparms);
    f_evals = 2*total_iters;
%
% TFQMR
%
elseif lmeth == 4 
    [step, errstep, total_iters] = dtfqmr(f0, f, x, gmparms);
    f_evals = 2*total_iters;
else
    error(' lmeth error in fdkrylov')
end
%
%
function z = dirder(x,w,f,f0)
% Finite difference directional derivative
% Approximate f'(x) w
% 
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

%
% Use a hardwired difference increment.
%
epsnew = 1.d-7;
%
n = length(x);
%
% scale the step
%
if norm(w) == 0
    z = zeros(n,1);
return
end
%
% Now scale the difference increment.
%
xs=(x'*w)/norm(w);
if xs ~= 0.d0
     epsnew=epsnew*max(abs(xs),1.d0)*sign(xs);
end
epsnew=epsnew/norm(w);
%
% del and f1 could share the same space if storage
% is more important than clarity.
%
del = x+epsnew*w;
f1 = feval(f,del);
z = (f1 - f0)/epsnew;
%
%
function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
% GMRES linear equation solver for use in Newton-GMRES solver
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any 
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%            params(3) (Optional) = reorthogonalization method
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarted GMRES
%              will be used as the linear solver.
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires givapp.m, dirder.m 

%
% initialization
%
errtol = params(1);
kmax = params(2);
reorth = 1;
if length(params) == 3
    reorth = params(3);
end
%
% The right side of the linear equation for the step is -f0. 
%
b = -f0;
n = length(b);
%
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
%
x = zeros(n,1);
r = b;
if nargin == 5
    x = xinit;
    r = -dirder(xc, x, f, f0)-f0;
end
%
%
h = zeros(kmax);
v = zeros(n,kmax);
c = zeros(kmax+1,1);
s = zeros(kmax+1,1);
rho = norm(r);
g = rho*eye(kmax+1,1);
errtol = errtol*norm(b);
error = [];
%
% Test for termination on entry.
%
error = [error,rho];
total_iters = 0;
if(rho < errtol) 
%   disp(' early termination ')
return
end
%
%
v(:,1) = r/rho;
beta = rho;
k = 0;
%
% GMRES iteration
%
while((rho > errtol) & (k < kmax))
    k = k+1;
%
%   Call directional derivative function.
%
    v(:,k+1) = dirder(xc, v(:,k), f, f0);
    normav = norm(v(:,k+1));
%
%   Modified Gram-Schmidt
%
    for j = 1:k
        h(j,k) = v(:,j)'*v(:,k+1);
        v(:,k+1) = v(:,k+1)-h(j,k)*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
    normav2 = h(k+1,k);
%
%   Reorthogonalize?
%
if  (reorth == 1 & normav + .001*normav2 == normav) | reorth ==  3
    for j = 1:k
        hr = v(:,j)'*v(:,k+1);
	h(j,k) = h(j,k)+hr;
        v(:,k+1) = v(:,k+1)-hr*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
end
%
%   Watch out for happy breakdown.
%
    if(h(k+1,k) ~= 0)
    v(:,k+1) = v(:,k+1)/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation.
%
    if k > 1
        h(1:k,k) = givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
%
%   Don't divide by zero if solution has  been found.
%
    nu = norm(h(k:k+1,k));
    if nu ~= 0
%        c(k) = h(k,k)/nu;
        c(k) = conj(h(k,k)/nu);
        s(k) = -h(k+1,k)/nu;
        h(k,k) = c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k) = 0;
        g(k:k+1) = givapp(c(k),s(k),g(k:k+1),1);
    end
%
%   Update the residual norm.
%
    rho = abs(g(k+1));
    error = [error,rho];
%
%   end of the main while loop
%
end
%
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
%
y = h(1:k,1:k)\g(1:k);
total_iters = k;
x = x + v(1:n,1:k)*y;
%
%
function vrot = givapp(c,s,vin,k)
%  Apply a sequence of k Givens rotations, used within gmres codes.
% 
%  C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
%  function vrot = givapp(c, s, vin, k)
%
vrot = vin;
for i = 1:k
    w1 = c(i)*vrot(i)-s(i)*vrot(i+1);
%
%   Here's a modest change that makes the code work in complex
%   arithmetic. Thanks to Howard Elman for this.
%
%    w2 = s(i)*vrot(i)+c(i)*vrot(i+1);
    w2 = s(i)*vrot(i)+conj(c(i))*vrot(i+1);
    vrot(i:i+1) = [w1,w2];
end
%
%
function [x, error, total_iters] = ...
                     dcgstab(f0, f, xc, params, xinit)
% Forward difference Bi-CGSTAB solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters]
%                    = dcgstab(f0, f, xc, params, xinit)
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarts are needed.
%
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires: dirder.m
%

%
% initialization
%
b = -f0; 
n = length(b); errtol = params(1)*norm(b); kmax = params(2); error = []; 
rho = zeros(kmax+1,1);
%
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
%
x = zeros(n,1);
r = b;
if nargin == 5
    x = xinit;
    r = -dirder(xc, x, f, f0)-f0;
end
%
hatr0 = r;
k = 0; rho(1) = 1; alpha = 1; omega = 1;
v = zeros(n,1); p = zeros(n,1); rho(2) = hatr0'*r;
zeta = norm(r); error = [error,zeta];
%
% Bi-CGSTAB iteration
%
while((zeta > errtol) & (k < kmax))
    k = k+1;
    if omega == 0
       error('Bi-CGSTAB breakdown, omega = 0');
    end
    beta = (rho(k+1)/rho(k))*(alpha/omega);
    p = r+beta*(p - omega*v);
    v = dirder(xc,p,f,f0);
    tau = hatr0'*v;
    if tau == 0
        error('Bi-CGSTAB breakdown, tau = 0');
    end 
    alpha = rho(k+1)/tau;
    s = r-alpha*v; 
    t = dirder(xc,s,f,f0);
    tau = t'*t;
    if tau == 0
       error('Bi-CGSTAB breakdown, t = 0');
    end
    omega = t'*s/tau; 
    rho(k+2) = -omega*(hatr0'*t);
    x = x+alpha*p+omega*s;
    r = s-omega*t;
    zeta = norm(r);
    total_iters = k;
    error = [error, zeta];
end
%
%

function [x, error, total_iters] = ...
                     dtfqmr(f0, f, xc, params, xinit)
% Forward difference TFQMR solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters]
%                    = dtfqmr(f0, f, xc, params, xinit)
%
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarts are needed.
%
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires: dirder.m
%

%
% initialization
%
b = -f0;
n = length(b); errtol = params(1)*norm(b); kmax = params(2); error = []; 
x = zeros(n,1);
r = b;
if nargin == 5
    x = xinit;
    r = -dirder(xc, x, f, f0)-f0;
end
%
u = zeros(n,2); y = zeros(n,2); w = r; y(:,1) = r; 
k = 0; d = zeros(n,1); 
v = dirder(xc, y(:,1),f,f0);
u(:,1) = v;
theta = 0; eta = 0; tau = norm(r); error = [error,tau];
rho = tau*tau;
%
% TFQMR iteration
%
while( k < kmax)
    k = k+1;
    sigma = r'*v;
%
    if sigma == 0
        error('TFQMR breakdown, sigma = 0')
    end
%
    alpha = rho/sigma;
%
% 
%
    for j = 1:2
%
%   Compute y2 and u2 only if you have to
%
        if j == 2 
            y(:,2) = y(:,1)-alpha*v;
            u(:,2) = dirder(xc, y(:,2),f,f0);
        end
        m = 2*k-2+j;
        w = w-alpha*u(:,j);
        d = y(:,j)+(theta*theta*eta/alpha)*d;
        theta = norm(w)/tau; c = 1/sqrt(1+theta*theta);
        tau = tau*theta*c; eta = c*c*alpha;
        x = x+eta*d;
%
%   Try to terminate the iteration at each pass through the loop
%
        if tau*sqrt(m+1) <=  errtol
            error = [error, tau];
            total_iters = k;
            return
        end
    end
%
%
%
    if rho == 0
        error('TFQMR breakdown, rho = 0')
    end
%
    rhon = r'*w; beta = rhon/rho; rho = rhon;
    y(:,1) = w + beta*y(:,2);
    u(:,1) = dirder(xc, y(:,1),f,f0);
    v = u(:,1)+beta*(u(:,2)+beta*v);
    error = [error, tau];
    total_iters = k;
end
%
%
