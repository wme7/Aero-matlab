% learing to use parfor
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 4
end
