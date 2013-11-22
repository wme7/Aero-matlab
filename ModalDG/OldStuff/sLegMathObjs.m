% Build Math objects for scaled Legendre polynomials. 

% M matrix
Mcoef = [1 1/12 1/180 1/2800 1/44100 1/698544 1/11099088 1/176679360];
M = diag(Mcoef(1:k+1));

% invM matrix
invM = inv(M);

% D matrix
Dcoef = [ ...
    0, 1, 0, (1/10), 0, (1/126), 0, (1/1716); ...
    0, 0, (1/6), 0, (1/70), 0, (1/924), 0; ...
    0, 0, 0, (1/60), 0, (1/756), 0, (1/10296); ...
    0, 0, 0, 0, (1/700), 0, (1/9240), 0; ...
    0, 0, 0, 0, 0, (1/8820), 0, (1/120120); ...
    0, 0, 0, 0, 0, 0, (1/116424), 0; ...
    0, 0, 0, 0, 0, 0, 0, (1/1585584); ...
    0, 0, 0, 0, 0, 0, 0, 0];
D = Dcoef(1:k+1,1:k+1);

% Scaled Legendre polynomials of deg 'l' evaluated at x = +1/2
Ln = zeros(k+1,1); % column vector
for l = 0:k; % Polynomials degree
    Ln(l+1) = sLegendreP(l,0.5);
end

% Scaled Legendre polynomials of deg 'l' evaluated at x = -1/2
Lp = zeros(k+1,1); % column vector
for l = 0:k; % Polynomials degree
    Lp(l+1) = sLegendreP(l,-0.5);
end