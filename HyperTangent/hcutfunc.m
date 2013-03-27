function h = hcutfunc(x,a,b)
% Compute Cut function scaling Hyperbolic tangent to the desired range.
% NOTE: tanh is constrained to work in the local range of [-4,4], we call
% this local domain $\xi$.
% coded by Manuel Diaz, NTU, IAM, 2013.03.27.

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
xi = (4-(-4))/(b-a)*(x_buffer - a) - 4;
h_tanh = 1/2*(tanh(xi)+1);

% define cut function h in all region
h(x_id_Euler) = 0;
h(x_id_Buffer) = h_tanh(:);
h(x_id_BGK) = 1;