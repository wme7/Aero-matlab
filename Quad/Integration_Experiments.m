clear;
% some integration experiments

% Example 1: E[p(x)], x~N(0,1), p(x) polynomial
degree = 10;
p = rand(1,degree+1); % make some polynomial
fprintf('adaptive quad, tol %.1g = %.10g\n', ...
        1e-10, quad(@(x) polyval(p,x).*normpdf(x), ...
                    -30,30,1e-10));

fprintf(' gauss-hermite quadrature\n');
for n=1:((degree+1)/2+4)
  % use hermite -- will be exact for n>=(degree+1)/2
  int = gaussHermite(n);
  % notice the change of variables!
  ep = polyval(p,int.x*sqrt(2))'*int.w/sqrt(pi);
  if (n==round((degree+1)/2))
    fprintf('--- the rest should be exact ----\n');
  end
  fprintf('n=%2d  E[p(x)] = %.10g\n',n,ep);
end

fprintf('\n monte carlo integration\n')
for n=1:6
  fprintf('n=10^%d  E[p(x)] = %.10g\n',n, ...
          mean(polyval(p,randn(10^n,1))));
end

pause();

% Example 2: E[p(x)|x>-0.1], x~N(0,1)
fprintf('adaptive quad, tol %.1g = %.10g\n', ...
        1e-10, quad(@(x) polyval(p,x).*normpdf(x), ...
                    -0.1,30,1e-10) / ...
        quad(@(x) normpdf(x),-0.1,30,1e-10));

fprintf('\n gauss-hermite quadrature (not the right rule)\n');
for n=1:2:60
  int = gaussHermite(n);
  m = int.x>-0.1;
  ep = polyval(p,int.x(m)*sqrt(2))'* ...
       int.w(m)/sqrt(pi) ...
       * sum(int.w)/sum(int.w(m));
  fprintf('n=%2d  E[p(x)] = %.10g\n',sum(m),ep);
end

fprintf('\n monte carlo integration\n')
for n=1:6
  x = randn(10^n,1);
  x = x(x>-0.1);
  fprintf('n=10^%d  E[p(x)] = %.10g\n',numel(x), ...
          mean(polyval(p,x)));
end
pause();

% Example 2: E[p(x)|x>-0.1], x~N(0,1)
% Example 2: E[p(x)|x>-0.1], x~N(0,1)
fprintf('adaptive quad, tol %.1g = %.10g\n', ...
        1e-10, quad(@(x) exp(x).*normpdf(x), ...
                    -30,30,1e-10));

fprintf('\n gauss-hermite quadrature (not the right rule)\n');
for n=1:20
  int = gaussHermite(n);
  ep = exp(int.x*sqrt(2))'*int.w/sqrt(pi);
  fprintf('n=%2d  E[p(x)] = %.10g\n',n,ep);
end

fprintf('\n monte carlo integration\n')
for n=1:6
  x = randn(10^n,1);
  fprintf('n=10^%d  E[p(x)] = %.10g\n',numel(x), ...
          mean(exp(x)));
end
pause();