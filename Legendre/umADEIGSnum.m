function  umADEIGSnum(p)

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

form = 'b*g+roc<m>y^ks';

hold off; clf;
for ntheta = 0:100
    
    % create symbolic theta variable
    theta = ntheta*2*pi/(101);
    
    % multiplying factor
    phase = exp(-sqrt(-1)*theta);
   
    % homework should have read (dx/2)*M\(D + F + G*exp(-i*theta))
    A = (1/2)*M\(D + F + G*phase);
    
    eA = eig(A);
    
    for j=1:p+1
        subplot(1,3,1); hold on; plot(theta,imag(eA(j)),   'r*'); hold off; title('real(\omega)');
        subplot(1,3,2); hold on; plot(theta,real(eA(j)), 'g+'); hold off;  title('imag(\omega)');
        subplot(1,3,3); hold on; plot(-real(eA(j)), -imag(eA(j)), 'bo'); title('-\omega'); hold off; 
    end
    
end
