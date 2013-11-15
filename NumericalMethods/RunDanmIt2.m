function 

Mass = rand(1,10000);
Length = rand(1,10000);
Width = rand(1,10000);
Height = rand(1,10000);

a=1;
tic
for I = 1:10000;
   Density(I) = Mass (I)/(Length(I)*Width(I)*Height(I));
end
toc

%% vectorization
tic
Density = Mass./(Length.*Width.*Height);
toc