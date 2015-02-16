%% Compute the sqrt of an Integer.

% 1. ask for integeter value
v = input('Enter a integer value: ');
%disp(v);

isInteger = 0;

while isInteger == 0;
    
    roundv = round(v);
    remain = roundv-v;
    if remain == 0
        disp('thx for following directions : )')
        isInteger = 1;
    else
        %disp('Please enter an integer');
        v = input('Enter a integer value (Please): ');
    end
    
end

% 2. compute sqrt

x0 = 10;    % guess value

% f(x) = x^2 - v
% f'(x) = 2*x

%f0 = x0^2-v;    % f(x0)
%df0 = 2*x0;     % f'(x0)

%x1 = x0 - f0/df0; % 1st approximation

%disp(x1);

%f1 = x1^2-v;    % f(x0)
%df1 = 2*x1;     % f'(x0)

%x2 = x1 - f1/df1; % 2nd

%disp(x2);

%f2 = x2^2-v;    % f(x0)
%df2 = 2*x2;     % f'(x0)

%x3 = x2 - f2/df2; % 2nd

%disp(x3);

%x1 = x0 - f(x0)/df(x0);
%x2 = x1 - f(x1)/df(x1);
%x3 = x2 - f(x2)/df(x2);

for i = 1:10
    xnew = x0 - f(x0)/df(x0);
    
    % for the next step
    x0 = xnew;
end

% 3. return the answer

disp(['the square root is: ',num2str(xnew)])