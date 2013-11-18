%% Finite Element Methods
%
%  We use the linear finite element method for solving the Poisson equation
%  as an example to illustrate the main ingredients of finite element
%  methods.
%
% Reference: 
%
% * <http://www.math.uci.edu/~chenlong/226/Ch2FEM.pdf Introduction to Finite Element Methods>
% * <http://www.math.uci.edu/~chenlong/226/Ch3FEMcode.pdf Progamming of Finite Element Methods>

%% Variational formulation
%
% The classic formulation of the Poisson equation reads as
%
% $$ - \Delta u = f  \hbox{ in }  \Omega, \qquad u  = g_D  \hbox{ on }
% \Gamma _D,  \qquad  \nabla u\cdot n = g_N  \hbox{ on } \Gamma _N, $$
%
% where $\partial \Omega = \Gamma _D\cup \Gamma _N$ and $\Gamma _D\cap \Gamma _N=\emptyset$. 
% We assume $\Gamma _D$ is closed and $\Gamma _N$ open.
% 
% Denoted by $H_{g_D}^1(\Omega)=\{v\in L^2(\Omega), \nabla v\in L^2(\Omega)
% \hbox{ and } v|_{\Gamma _D} = g_D\}$. Multiplying the Poisson equation by
% a test function $v$ and using integration by parts, we obtain the weak
% formulation of the Poisson equation: find $u\in H_{g_D}^1(\Omega)$ such
% that for all $v\in H_{0_D}^1$:
%
% $$ a(u,v) := \int _{\Omega} \nabla u\cdot \nabla v\, {\rm dxdy} = \int _{\Omega} fv \, {\rm dxdy} + \int _{\Gamma _N} g_N v \,{dS}.$$
%
% Let $\mathcal T$ be a triangulation of $\Omega$. We define the linear
% finite element space on $\mathcal T$ as 
% 
% $$
% \mathcal V_{\mathcal T} = \{v\in C(\bar \Omega) : v|_{\tau}\in \mathcal P_k, \forall \tau \in \mathcal T\}. 
% $$
%
% where $\mathcal P_k$ is the polynomial space with degree $\leq k$. 
%
% The finite element method for solving the Poisson
% equation is to find $u\in \mathcal V_{\mathcal T}\cap H_{g_D}^1(\Omega)$ 
% such that for all $v\in \mathcal V_{\mathcal T}\cap H_{0_D}^1(\Omega)$:
%
% $$ a(u,v) = \int _{\Omega} fv \, {\rm dxdy} + \int _{\Gamma _N} g_N v \,{dS}.$$
%
%
%% Finite element space
%
% We take linear finite element spaces as an example. For each vertex $v_i$
% of $\mathcal T$, let $\phi _i$ be the piecewise linear function such that
% $\phi _i(v_i)=1$ and $\phi _i(v_j)=0$ when $j\neq i$. The basis function
% in 1-D and 2-D is illustrated below. It is also called hat function named
% after the shape of its graph.
% 
clf;
x = 0:1/5:1;
u = zeros(length(x),1);
u(2) = 1;
figure(1)
set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.4,0.3]);
subplot(1,2,1); hold on; 
plot(x,0,'k.','MarkerSize',18); 
plot(x,u,'-','linewidth',1.2);
subplot(1,2,2); hold on;
for k = 1:length(x)
    u = zeros(length(x),1); u(k) = 1;
    plot(x,0,'k.','MarkerSize',18); 
    plot(x,u,'-','linewidth',1.2);
end
%%
% 2-D hat basis
clf; set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.4,0.4]);
node = [0,0; 1,0; 1,1; 0,1];
elem = [2 3 1; 4 1 3];
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
u = zeros(size(node,1),1);
u(6) = 1;
showmesh(node,elem,'facecolor','none'); hold on;
showsolution(node,elem,u,[30,26],'facecolor','g','facealpha',0.5,'edgecolor','k');
%%
% Then it is easy to see $\mathcal V_{\mathcal T}$ is spanned by $\{\phi
% _i\}_{i=1}^{N}$ and thus for a finite element function $v=\sum
% _{i=1}^Nv_i\phi _i$.  


%% Matrix assembeling


%% Right hand side

%% Boundary condition
