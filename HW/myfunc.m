function ME = myfunc(x)
% Break Test

for i = 1:100
    try
        x = x - 1;
        if x < 0
            error('x has become negative')
        end
    catch ME
        fprintf('x become negative')
        break
    end
end
    
    