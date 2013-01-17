function  A = umADEIGSwhole(Ncell,p)

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

% build entire matrix
Ntotal = Ncell*(p+1);
A    = zeros(Ntotal);

for cell = 1:Ncell
    id1 = ((cell-1)*(p+1)+1):cell*(p+1);
    A(id1,id1) = (1/2)*M\(D+F);
    
    if(cell~=1)
        id2 = id1-(p+1);
    else
        id2 = id1+(Ncell-1)*(p+1);
    end
    
    A(id1, id2) = (1/2)*M\G;
end 

eA = eig(A);

plot(-real(eA), -imag(eA), 'bo'); 

