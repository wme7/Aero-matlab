
for p=0:140
    x = sin(0.5*pi*linspace(-1,1,p+1))';
    
    V = LEGvdm(x, p);
    
    cgraph(p+1) = cond(V);
    
end

plot(cgraph, 'r*'); 
hold on; plot(cgraph, 'k-');  hold off;
xlabel('p-order'); ylabel('cond(V)');
    