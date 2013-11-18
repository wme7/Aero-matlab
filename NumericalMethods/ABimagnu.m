
% compute AB2 roots on imaginary axis i*[0,1]
N = 100;
alpha = linspace(0, 1, N);
for n=1:N
    nu = i*alpha(n);
    p = [1,-1-3*nu/2, nu/2];
    ro(n) = max(abs(roots(p)));
end
plot(alpha, sign(ro-1).*abs(log10(abs(ro-1))), 'r-')
 
% compute AB3 roots on imaginary axis i*[0,1]

for n=1:N
    nu = i*alpha(n);
    p = [1,-1-(23*nu/12), (16*nu/12), -(5*nu/12)];
    ro(n) = max(abs(roots(p)));
end
hold on;
plot(alpha, sign(ro-1).*abs(log10(abs(ro-1))), 'b-')
hold off;

legend('AB2', 'AB3');

xlabel('|nu|');
ylabel('sign(root-1)*abs(log10(|maximum(|root|)-1|))', 'FontSize', 16);