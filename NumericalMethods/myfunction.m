function f = myfunction(x)
%% Function for Exercise problem:
% (1.) Just for this case if x = 0 then must be f(x)=0,
% (2.) m = is predefinet to evaluate the cases m={2,9,5.5}.
m=9;
if x==0
    f=0;
else
    f=(exp(-(1/x-1))*(1-x)^(m-1))/(x^(m+1));
end
