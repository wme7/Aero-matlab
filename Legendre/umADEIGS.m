function  A = umADEIGS(p)

% create symbolic theta variable
theta = sym('theta');

% multiplying factor
phase = exp(-sqrt(-1)*theta);

Dhat = zeros(p+1);
  for n=0:p
    for m=n+1:2:p
      Dhat(n+1,m+1) =  (2*n+1);
    end
  end

  % M = int_{-1}^{1} L_n(x)*L_m(x) dx (note dx=1 assumed)
  M = diag(2./(2*(0:p)+1));

  % D_nm = int_{-1}^{1}  L_n(x)*L_m(x)*dx
  D = M*Dhat;

   % F_nm = L_n(-1)*L_m(-1)
  F = (-1).^[0:p]'*(-1).^[0:p];

  % G_nm = -L_n(-1)*L_m(1)
  G = -(-1).^[0:p]'*ones(1,p+1);

  % homework should have read (dx/2)*M\(D + F + G*exp(-i*theta))
  A = (1/2)*M\(D + F + G*phase);
  
 
  eA = eig(A)
  eA = simplify(eA)

  % evaluate the Taylor series for each eigenvalue
  for j=1:p+1
      
      mode = j
   
      ev = eA(j);
      theta = 0;
      
      for k=1:10
          coeff(k) = eval(ev)/(factorial(k-1));
          ev = diff(ev);
      end
      
      coeff'
  end
  
  eV = zeros(p+1, 100);
  for j=0:99
      theta = j*2*pi/100;
      
      eV(:,j+1) = eval(eA);
  end
  
  theta = linspace(0,2*pi,100);
  form = 'bgrcmyk';
  hold off; clf;
  for j=1:p+1
      
    subplot(1,2,1); hold on;      plot(theta, imag(eV(j,:)), form(j)); hold off;  title('Im( i*\omega*dx)')
    subplot(1,2,2); hold on;      plot(theta, real (eV(j,:)), form(j)); hold off;   title('Re( i*\omega*dx)')
  end
  
