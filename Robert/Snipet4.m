% Snipet 4

fprintf('Fibonacci Series \n')

fmm=1; % f_{n-2} 
fm =1; % f_{n-1}
fprintf('number: %f \n',fmm)
fprintf('number: %f \n',fm )
for n = 3:10
    
    fn = fm + fmm;
    fprintf('number: %f \n',fn)
    
    % preparing for the next step
    fmm = fm; % f_{n-2} <- f_{n-1}
    fm  = fn; % f_{n-1} <- f_{n}
    
end
fprintf('\n')