% 2D Lattice Boltzmann (BGK) model of a fluid.
%  c4  c3   c2    D2Q9 model. At each timestep, particle densities propagate
%    \  |  /      outwards in the directions indicated in the figure. An
%  c5 -c9 - c1    equivalent 'equilibrium' density is found, and the densities
%    /  |  \      relax towards that state, in a proportion governed by omega.
%  c6  c7   c8
%                  Iain Haslam, April 2006.
omega=1.0; density=1.0; t1=4/9; t2=1/9; t3=1/36; c_squ=1/3; nx=1; ny=101;
F=repmat(density/9,[nx ny 9]); F_EQ=F; msize=nx*ny; qi=[0:msize:msize*7];
bound=zeros(nx,ny); bound(:,[1 ny])=1;%bound([14:19],[14:16])=1;
obs=find(bound); 
bounds_to_bounce=[obs+qi(1) obs+qi(2) obs+qi(3) obs+qi(4) ... 
	obs+qi(5) obs+qi(6) obs+qi(7) obs+qi(8)];
bounds_bounced=[obs+qi(5) obs+qi(6) obs+qi(7) obs+qi(8) ...
	obs+qi(1) obs+qi(2) obs+qi(3) obs+qi(4)];
avu=1; prevavu=1; ts=0; deltaU=1e-9; numactivenodes=sum(sum(1-bound));
while (ts<100000 & 1e-10<abs((prevavu-avu)/avu)) | ts<100
    ts=ts+1;
	% Propagate
	F(:,:,4)=F([2:nx 1],[ny 1:ny-1],4);F(:,:,3)=F(:,[ny 1:ny-1],3);
	F(:,:,2)=F([nx 1:nx-1],[ny 1:ny-1],2);F(:,:,5)=F([2:nx 1],:,5);
	F(:,:,1)=F([nx 1:nx-1],:,1);F(:,:,6)=F([2:nx 1],[2:ny 1],6);
	F(:,:,7)=F(:,[2:ny 1],7); F(:,:,8)=F([nx 1:nx-1],[2:ny 1],8);
	BOUNCEDBACK=F(bounds_to_bounce); %Densities bouncing back at next timestep
	% Relax; calculate equilibrium state with equivalent speed and density to F 
	DENSITY = sum(F,3);	
	UX=(sum(F(:,:,[1 2 8]),3)-sum(F(:,:,[4 5 6]),3))./DENSITY;
	UY=(sum(F(:,:,[2 3 4]),3)-sum(F(:,:,[6 7 8]),3))./DENSITY;
	%UX(1,2:ny-1)=UX(1,2:ny-1)+deltaU;
	UX=UX+deltaU;
	UX(obs)=0; UY(obs)=0; DENSITY(obs)=0;
	U_SQU=UX.^2+UY.^2; U_C2=UX+UY; U_C4=-UX+UY; U_C6=-U_C2; U_C8=-U_C4;
	% stationary
	F_EQ(:,:,9)=t1*DENSITY.*(1-U_SQU/(2*c_squ));
	% nearest-neighbours
	F_EQ(:,:,1)=t2*DENSITY.*(1+UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,3)=t2*DENSITY.*(1+UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,5)=t2*DENSITY.*(1-UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,7)=t2*DENSITY.*(1-UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
	% next-nearest neighbours
	F_EQ(:,:,2)=t3*DENSITY.*(1+U_C2/c_squ+0.5*(U_C2/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,4)=t3*DENSITY.*(1+U_C4/c_squ+0.5*(U_C4/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,6)=t3*DENSITY.*(1+U_C6/c_squ+0.5*(U_C6/c_squ).^2-U_SQU/(2*c_squ));
	F_EQ(:,:,8)=t3*DENSITY.*(1+U_C8/c_squ+0.5*(U_C8/c_squ).^2-U_SQU/(2*c_squ));
	F=omega*F_EQ+(1-omega)*F;
	F(bounds_bounced)=BOUNCEDBACK;
	prevavu=avu;avu=sum(sum(UX))/numactivenodes; ts=ts+1;
end
%The bounceback boundary condition means that the actual boundary lies
%mid-way between the open and closed nodes. The width of the channel is
%therefore ny-2.
L=(ny-2)/2; arrindices=[-(L-0.5):(L-0.5)]; EXACT=-(arrindices.^2-L^2);
tau = 1/omega;
%Calculate the analytical solution's peak velocity
temp_F = density * deltaU / tau;
temp_viscosity = (tau-0.5)/3;
calcfactor=temp_F/(2*temp_viscosity);
CALC=[UX(1,2:ny-1)];
EXACT=EXACT*calcfactor;
figure;plot(arrindices,abs(1-(EXACT./CALC)));
maxerror=max(abs(1-(EXACT./CALC)));
title('Relative error of LBM solution across channel width');
fprintf('RESULTS\n')
fprintf('ts=%d, maxerror=%g\n',ts, maxerror);
%Calculate relative proportional error in centre of channel
[maxcalc, maxindex] = max(CALC);
maxerrorcentre = (EXACT(maxindex)-CALC(maxindex))/EXACT(maxindex)
figure;plot(arrindices,EXACT);hold on; plot(arrindices,CALC,'o');
title('Comparison of analytical with LBM results');
xlabel('location in channel'); ylabel('speed')
