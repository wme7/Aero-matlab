% r+1'th order, r+1 stage RK reproduces the first
% r+1 terms of the Taylor expansion in time
% 
% y^{n+1} = {\sum_{m=0}^{m=r} ((mu)^r)/(r!) }*y^n 
% 

function RKabstab(p, form)

  % set up polynomial coefficients 
  % (highest order coefficient first)
  mypoly = 1./factorial(p:-1:0)';
  
  N = 100;

  rootstore = [];

  % plot dependence of polynomial roots on theta at N values
  for ntheta = 1:N
    theta = 2*pi*(ntheta-1)/(N-1);
    
    mypoly(p+1) = 1 - exp(i*theta);
    
    ro      = roots(mypoly);
    
    rootstore = [rootstore; ro];
  end
  
  % sort roots by argument to plot curve
  angles = atan2(imag(rootstore), real(rootstore)+1);
  [dummy,order] = sort(angles);
  rootstore = rootstore(order); 

  % plot sorted roots
  hold on; 
  plot(real(rootstore), imag(rootstore), form); 

  % plot real and imaginary axes
  plot([-5 5],[ 0 0], 'k-');
  plot([ 0 0],[-5 5], 'k-');
  axis([-3 0.5 -3.25 3.25]); 
  xlabel('real(alpha)'); ylabel('imag(alpha)');

