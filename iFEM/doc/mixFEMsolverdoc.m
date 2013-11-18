%% Mixed FEM solver

%% Uzawa PCG
parabolicMixedFEM('h^2','uzawapcg');
parabolicMixedFEM('1','uzawapcg');

%% Transfer MG
parabolicMixedFEM('h^2','dmg');
parabolicMixedFEM('1','dmg');