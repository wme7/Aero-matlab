for p=0:140
    x = sin(0.5*pi*linspace(-1,1,p+1))'; % Compute Chebishev nodes
    
    V = legendreVDM(x,p);
    
    cgraph(p+1) = cond(V); % condition matlab test using L2-norm.
    
end

plot(cgraph,'r*'); 
hold on; plot(cgraph,'k-');  hold off;
xlabel('p-order'); ylabel('cond(V)');
    