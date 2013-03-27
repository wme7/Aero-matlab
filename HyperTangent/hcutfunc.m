function h = hcutfunc(x,a,b)
% Compute Cut function scaling Hyperbolic tangent to the desired range.
% NOTE: 'tanh' is constrained to work in the local range of [-4,4] and
% 'cos' is constrained to work in the local range of [pi,2*pi]. We refer 
% to this local domain $\xi$.
% coded by Manuel Diaz, NTU, IAM, 2013.03.27.
modelfunc = 2;

% Idenfity regions,
h_Euler = (x<a);
h_BGK = (x>=b);
h_Buffer = (x>=a & x<b);

% Identify node Id for every region,
x_id_Euler = find(h_Euler);
x_id_BGK = find(h_BGK);
x_id_Buffer = find(h_Buffer);

% define h in buffer region
x_buffer = x(x_id_Buffer);
a = x_buffer(1); b = x_buffer(end);

% apply transformation
switch modelfunc
    case{1} % using tanh(x) for xi = [-4,4]
        xi = (4-(-4))/(b-a)*(x_buffer - a) - 4;
        h_func = 1/2*(tanh(xi)+1);
    case{2} % using cos(x) for xi = [pi,2*pi]
        xi = (2*pi-(pi))/(b-a)*(x_buffer - a) - pi;
        h_func = 1/2*(cos(xi)+1);
end
% define cut function h in all region
h(x_id_Euler) = 0;
h(x_id_Buffer) = h_func(:);
h(x_id_BGK) = 1;