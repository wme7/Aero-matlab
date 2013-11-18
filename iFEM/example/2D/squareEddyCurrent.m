pde = eddycurrentdata1;
[node,elem] = squaremesh([0,1,0,1],1/16);
bdFlag = setboundary(node,elem,'Dirichlet');
err = zeros(4,1); N = zeros(4,1);
for i=1:4
    [node,elem,bdFlag]=uniformrefine(node,elem,bdFlag);
    option.solver = 'cg';
    [u,edge,eqn] = eddycurrent(node,elem,pde,bdFlag,option);
    err(i) = getHdiverrorRT0(node,elem,0,u);
    N(i) = size(u,1);
end
showrate(N,err,2);