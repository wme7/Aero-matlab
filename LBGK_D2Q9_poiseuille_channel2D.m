% Gianni Schena  July 2005, schena@units.it
% Lattice Boltzmann LBE, geometry: D2Q9, model: BGK
% Application to permeability in porous media 

Restart=false % to restart from an earlier convergence
logical(Restart);

if Restart==false;
close all, clear all % start from scratch and clean ...
Restart=false;
% type of channel geometry ; 
% one of the flollowing flags == true
Pois_test=true, % no obstacles in the 2D channel
% porous systems
obs_regolare=false % 
obs_irregolare=false % 
tic
%   IN
% |vvvv|    + y
% |vvvv|     ^
% |vvvv|     | -> + x
%  OUT

% Pores in 2D : Wet and Dry locations (Wet ==1 , Dry ==0 )
wXh_Dry=[3,1];wXh_Wet=[3,4];

if obs_regolare, % with internal obstacles 
    
A=repmat([zeros(wXh_Dry),ones(wXh_Wet)],[1,3]);A=[A,zeros(wXh_Dry)];
B=ones(size(A)); 
C=[A;B]  ; D=repmat(C,4,1);
D=[B;D]
end

if obs_irregolare, % with int obstacles 
A1=repmat([zeros(wXh_Dry),ones(wXh_Wet)],[1,3]); 
A1=[A1,zeros(wXh_Dry)]  ;
B=ones(size(A1)); 
C1=repmat([ones(wXh_Wet),zeros(wXh_Dry)],[1,3]); C1=[C1,ones(wXh_Dry)];
E=[A1;B;C1;B]; 
D=repmat(E,2,1);
D=[B;D]
end

if ~Pois_test
figure,imshow(D,[]) 
Channel2D=D;
Len_Channel_2D=size(Channel2D,1); % Length
Width=size(Channel2D,2); % should not be hod
Channel_2D_half_Width=Width/2,
end

% test without obstacles (i.e. 2D channel & no obstacles)

if Pois_test
%over-writes the definition of the pore space
clear Channel2D
Len_Channel_2D=36, % lunghezza canale 2d
Channel_2D_half_Width=8; Width=Channel_2D_half_Width*2;
Channel2D=ones(Len_Channel_2D,Width); % define wet area
%Channel2D(6:12,6:8)=0; % put fluid obstacle
imshow(Channel2D,[]);
end

[Nr Mc]=size(Channel2D); % Number rows and Munber columns

% porosity
porosity=nnz(Channel2D==1)/(Nr*Mc)


% FLUID PROPERTIES
% physical properties
cs2=1/3; % 
cP_visco=0.5; % [cP] 1 CP Dinamic water viscosity 20 C
density=1.; % fluid density 
Lky_visco=cP_visco/density; % lattice kinematic viscosity 
omega=(Lky_visco/cs2+0.5).^-1; %  omega: relaxation frequency
%Lky_visco=cs2*(1/omega - 0.5) , % lattice kinematic viscosity
%dPdL= Pressure / dL;% External pressure gradient [atm/cm]

uy_fin_max=-0.2; 
%dPdL = abs( 2*Lky_visco*uy_fin_max/(Channel_2D_half_Width.^2) ); 
dPdL=-0.0125;
uy_fin_max=dPdL*(Channel_2D_half_Width.^2)/(2*Lky_visco); % Poiseuille Gradient;
% max poiseuille final  velocity on the flow profile
uy0=-0.001; ux0=0.0001; %  linear vel .. inizialization

% 
% uy_fin_max=-0.2; % max poiseuille final  velocity on the flow profile
% omega=0.5, cs2=1/3; % omega: relaxation frequency
% Lky_visco=cs2*(1/omega - 0.5) , % lattice kinematic viscosity
% dPdL = abs( 2*Lky_visco*uy_fin_max/(Channel_2D_half_Width.^2) ); % Poiseuille Gradient;
% 

uyf_av=uy_fin_max*(2/3);; % average fluid velocity on the profile

x_profile=([-Channel_2D_half_Width:+Channel_2D_half_Width-1]+0.5);
uy_analy_profile=uy_fin_max.*(1-  ( x_profile /Channel_2D_half_Width).^2 ); % analytical velocity profile

av_vel_t=1.e+10; % inizialization (t=0)
%PixelSize= 5; % [Microns]
%dL=(Nr*PixelSize*1.0E-4); % sample hight [cm]


%
% EXPERIMENTAL SET-UP
% inlet and outlet buffers
inb=2, oub=2; % inlet and outlet buffers thickness
% add fluid at the inlet (top) and outlet (down)
inlet=ones(inb,Mc); outlet=ones(oub,Mc);
Channel2D=[ [inlet]; Channel2D ;[outlet] ] ; % add flux in and down (E to W)
[Nr Mc]=size(Channel2D); % update size
% boundaries related to the experimental set up
wb=2; % wall thickness
Channel2D=[zeros(Nr,wb), Channel2D , zeros(Nr,wb)]; % add walls (no fluid leak)
[Nr Mc]=size(Channel2D); % update size
uy_analy_profile=[zeros(1,wb), uy_analy_profile, zeros(1,wb) ] ; % take into account walls
x_pro_fig=[[x_profile(1)-[wb:-1:1]], [x_profile, [1:wb]+x_profile(end)] ];

% Figure plots analytical parabolic profile
figure(20), plot(x_pro_fig,uy_analy_profile,'-'), grid on,
title('Analytical parab. profile for Poiseuille planar flow in a channel')


% VISUALIZE PORE SPACE & FLUID OSTACLES & MEDIAL AXIS
figure, imshow(Channel2D); title('Vassel geometry');
Channel2D=logical(Channel2D);
% obstacles for Bounce Back ( in front of the grain)
Obstacles=bwperim(Channel2D,8); % perimeter of the grains for bounce back Bound.Cond.
border=logical(ones(Nr,Mc));
border([1:inb,Nr-oub:Nr],[wb+2:Mc-wb-1])=0;
Obstacles=Obstacles.*(border);
figure, imshow(Obstacles); title(' Fluid obstacles (in the fluid)' );
% 
Medial_axis=bwmorph(Channel2D,'thin',Inf); %
figure, imshow(Medial_axis); title('Medial axis');
figure(10) % used to visualize evolution of rho
figure(11) % used to visualize ux
figure(12) % used to visualize uy (i.e. top -> down)

% INDICES
% Wet locations etc.
[iabw1 jabw1]=find(Channel2D==1); % indices i,j, of active lattice locations i.e. pore
lena=length(iabw1); % number of active location i.e. of pore space lattice cells
ija= (jabw1-1)*Nr+iabw1; % equivalent single index (i,j)->> ija for active locations
% absolute (single index) position of the obstacles in for bounce back in Channel2D
% Obstacles 
[iobs jobs]=find(Obstacles);lenobs=length(iobs); ijobs= (jobs-1)*Nr+iobs; % as above
% Medial axis of the pore space
[ima jma]=find(Medial_axis); lenma=length(ima);  ijma= (jma-1)*Nr+ima; % as above
% Internal wet locations : wet & ~obstables
% (i.e. internal wet lattice location non in contact with dray locations)
[iawint jawint]=find(( Channel2D==1 & ~Obstacles)); % indices i,j, of active lattice locations
lenwint=length(iawint); % number of internal (i.e. not border) wet locations
ijaint= (jawint-1)*Nr+iawint; % equivalent singl
NxM=Nr*Mc;

% DIRECTIONS: E N W S NE NW SW SE ZERO (ZERO:Rest Particle)
%    y^
%  6 2 5           ^         NW  N  NE
%  3 9 1 ... +x-> +y         W   RP  E
%  7 4 8                     SW  S  SE
%   -y
% x & y components of velocities , +x is to est , +y is to nord
East=1; North=2; West=3; South=4; NE=5; NW=6; SW=7; SE=8; RP=9;
N_c=9 ; % number of directions
% versors D2Q9
C_x=[1 0 -1  0 1 -1 -1  1 0]; 
C_y=[0 1  0 -1 1  1 -1 -1 0]; C=[C_x;C_y]

% BOUNCE BACK SCHEME
% after collision the fluid elements densities f are sent back to the
% lattice node they come from with opposite direction
% indices opposite to 1:8 for fast inversion after bounce
ic_op = [3 4 1 2 7 8 5 6]; %   i.e. 4 is opposite to 2 etc.

% PERIODIC BOUNDARY CONDITIONS - reinjection rules
yi2=[Nr , 1:Nr , 1]; % this definition allows implemening Period Bound Cond
%yi2=[1, Nr , 2:Nr-1 , 1,Nr]; % re-inj the second last to as first
% directional weights (density weights)
w0=16/36. ; w1=4/36. ; w2=1/36.;
W=[ w1 w1 w1 w1 w2 w2 w2 w2 w0];
%c constants (sound speed related)
cs2=1/3; cs2x2=2*cs2; cs4x2=2*cs2.^2;
f1=1/cs2; f2=1/cs2x2; f3=1/cs4x2;
f1=3., f2=4.5; f3=1.5; % coef. of the f equil.

% declarative statemets
f=zeros(Nr,Mc,N_c); % array of fluid density distribution
feq=zeros(Nr,Mc,N_c); % f at equilibrium
rho=ones(Nr,Mc); % macro-scopic density
temp1=zeros(Nr,Mc);
ux=zeros(Nr,Mc);   uy=zeros(Nr,Mc); uyout=zeros(Nr,Mc);  % dimensionless velocities
uxsq=zeros(Nr,Mc); uysq=zeros(Nr,Mc);   usq=zeros(Nr,Mc);  % higher degree velocities

% initialization arrays : start values in the wet area
for ia=1:lena % stat values in the active cells only ; 0 outside
    i=iabw1(ia);  j=jabw1(ia);
    f(i,j,:)=1/9; % uniform density distribution for a start
end
uy(ija)=uy0; ux(ija)=ux0; % initialize fluid velocities
rho(ija)=density;

% EXTERNAL (Body) FORCES e.g. inlet pressure or inlet-outlet gradient
% directions: E N W S NE NW SW SE ZERO
force = -dPdL*(1/6)*1*[0 -1 0 1 -1 -1 1  1  0]'; %;
%...                   E  N E S NE NW SW SE RP ...
% the pressure pushes the fluid down i.e. N to S

% While .. MAIN TIME EVOLUTION LOOP
StopFlag=false; % i.e. logical(0)
Max_Iter=3000; % max allowed number of iteration
Check_Iter=1; Output_Every=20; % frequency of check & output
Cur_Iter=0; % current iteration counter inizialization
toler=1.0e-8; % tollerance to declare convegence
Cond_path=[]; % recording values of the convergence criterium
density_path=[]; % recording aver. density values for convergence
end % ends if restart

if(Restart==true)
 StopFlag=false;  Max_Iter=Max_Iter+3000; toler=1.0e-12; 
end


while(~StopFlag)
    Cur_Iter=Cur_Iter+1 % iteration counter update

    % density and moments
    rho=sum(f,3); % density

    if Cur_Iter >1 % use inizialization ux uy to start
        % Moments ... Note:C_x(9)=C_y(9)=0
        ux=zeros(Nr,Mc); uy=zeros(Nr,Mc);
        for ic=1:N_c-1;
            ux = ux + C_x(ic).*f(:,:,ic) ; uy = uy + C_y(ic).*f(:,:,ic)  ;
        end
       % uy=f(:,:,2) +f(:,:,5)+f(:,:,6)-f(:,:,4)-f(:,:,7)-f(:,:,8); % in short !
       % ux=f(:,:,1) +f(:,:,5)+f(:,:,8)-f(:,:,3)-f(:,:,6)-f(:,:,7); % in short !
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ux(ija)=ux(ija)./rho(ija); uy(ija)=uy(ija)./rho(ija);
    uxsq(ija)=ux(ija).^2; uysq(ija)=uy(ija).^2; 
    usq(ija)=uxsq(ija)+uysq(ija); %

    % weighted densities : rest particle, principal axis, diagonals
    rt0 = w0.*rho; rt1 = w1.*rho; rt2 = w2.*rho;
    
    % Equilibrium distribution
    % main  directions ( + cross)
    feq(ija)= rt1(ija) .*(1 +f1*ux(ija) +f2*uxsq(ija) -f3*usq(ija));
    feq(ija+NxM*(2-1))= rt1(ija) .*(1 +f1*uy(ija) +f2*uysq(ija) -f3*usq(ija));
    feq(ija+NxM*(3-1))= rt1(ija) .*(1 -f1*ux(ija) +f2*uxsq(ija) -f3*usq(ija));
    %feq(ija+NxM*(3)=f(ija)-2*rt1(ija)*f1.*ux(ija); % much faster... !!
    feq(ija+NxM*(4-1))= rt1(ija) .*(1 -f1*uy(ija) +f2*uysq(ija) -f3*usq(ija));
    
    % diagonals (X diagonals) (ic-1)
    feq(ija+NxM*(5-1))= rt2(ija) .*(1 +f1*(+ux(ija)+uy(ija)) +f2*(+ux(ija)+uy(ija)).^2 -f3.*usq(ija));
    feq(ija+NxM*(6-1))= rt2(ija) .*(1 +f1*(-ux(ija)+uy(ija)) +f2*(-ux(ija)+uy(ija)).^2 -f3.*usq(ija));
    feq(ija+NxM*(7-1))= rt2(ija) .*(1 +f1*(-ux(ija)-uy(ija)) +f2*(-ux(ija)-uy(ija)).^2 -f3.*usq(ija));
    feq(ija+NxM*(8-1))= rt2(ija) .*(1 +f1*(+ux(ija)-uy(ija)) +f2*(+ux(ija)-uy(ija)).^2 -f3.*usq(ija));
    % rest particle (.) ic=9
    feq(ija+NxM*(9-1))= rt0(ija) .*(1 - f3*usq(ija));

    %Collision (between fluid elements)omega=relaxation frequency
    f=(1.-omega).*f + omega.*feq;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %add external body force due to the pressure gradient prop. to dPdL
    for ic=1:N_c;%-1
        for ia=1:lena
            i=iabw1(ia);  j=jabw1(ia);
            % if Obstacles(i,j)==0 % the i,j is not aderent to the boundaries
            % if ( f(i,j,ic) + force(ic) ) >0; %! avoid negative distributions
            %i=1 ;% force only on the first row !
            f(i,j,ic)= f(i,j,ic) + force(ic);
            % end
            % end
        end
    end

   

    % % STREAM
    % Forward Propagation step & % Bounce Back (collision fluid with obstacles)
    %f(:,:,9) = f(:,:,9); % Rest element do not move
   
    feq = f; % temp storage of f in feq
        for ic=1:1:N_c-1, % select velocity layer

        ic2=ic_op(ic); % selects the layer of the velocity opposite to ic for BB
        temp1=feq(:,:,ic); %

        % from wet location that are NOT on the border to other wet locations
        for ia=1:1:lenwint % number of internal (i.e. not border) wet locations
            i=iawint(ia);  j=jawint(ia);  % so that we care for the wet space only !
            i2 = i+C_y(ic); j2 = j+C_x(ic); % Expected final locations to move
            i2=yi2(i2+1); % i2 corrected for PBC when necessary (flow out re-fed to inlet)
            % i.e the new position (i2,j2) is sure another wet location
            % therefore normal propagation from (i,j) to (i2,j2) on layer ic
            f(i2,j2,ic)=temp1(i,j); % see circshift(..) fnct for circularly shifts
        end ; % i and j single loop


        % from wet locations that ARE on the border of obstacles
        for ia=1:1:lenobs % wet border locations
            i=iobs(ia);  j=jobs(ia);  % so that we care for the wet space only !
            i2 = i+C_y(ic); j2 = j+C_x(ic); % Expected final locations to move
            i2=yi2(i2+1); % i2 corrected for PBC

            if( Channel2D(i2,j2) ==0 ) % i.e the new position (i2,j2) is dry
                f(i,j,ic2) =temp1(i,j); % invert direction: bounce-back in the opposite direction ic2
            else % otherwise, normal propagation from (i,j) to (i2,j2) on layer ic
                f(i2,j2,ic)=temp1(i,j); % see circshift(..) fnct for circularly shifts
            end ; % b.b. and propagations

        end ; % i and j single loop
        % special treatment for Corners
        %   f(1,wb+1,ic)=temp1(Nr,Mc-wb);      f(1,Mc-wb,ic)=temp1(Nr,wb+1);
        %   f(Nr,wb+1,ic)=temp1(1,Mc-wb);      f(Nr,Mc-wb,ic)=temp1(1,wb+1);

    end ; %  for ic direction

    % ends of Forward Propagation step &  Bounce Back Sections

    % re-calculate  uy as uyout for convergence
    rho=sum(f,3); % density
    % check velocity
    uyout= zeros(Nr,Mc);
    for ic=1:N_c-1;
        uyout= uyout + C_y(ic).*f(:,:,ic) ; % flow dim.less velocity out
    end
   % uyout(ija)=uyout(ija)./rho(ija); % from momentum to velocity

    % Convergence check on velocity values
    if (mod(Cur_Iter,Check_Iter)==0) ; % check for convergence every 'Check_Iter' iterations

        % variables monitored
        % mean density and
        vect=rho(ija); vect=vect(:); 
        cur_density=mean(vect);
        % mean 'interstitial' velocity
        % uy(ija)=uy(ija)/rho(ija); ?
        vect=uy(ija); av_vel_int= mean(vect)  ; % seepage velocity (in the wet area)
        % on the whole cross-sectional area of flow (wet + dry)
        av_vel_int=av_vel_int*porosity, % av. vel. on the wet + dry area
        %av_vel_int=mean2(uy),
        av_vel_tp1 = av_vel_int; 
        Condition=abs( abs(av_vel_t/av_vel_tp1 )-1), % should --> 0

        Cond_path=[Cond_path, Condition]; % records the convergence path (value)
        density_path=[density_path, cur_density];
        %
        av_vel_t=av_vel_tp1; % time t & t+1 

        if (Condition < toler) | (Cur_Iter > Max_Iter)
            StopFlag=true;
            display( 'Stop iteration: Convergence met or iteration exeeding the max allowed' )
            display( ['Current iteration: ',num2str(Cur_Iter),...
                ' Max Number of iter: ',num2str(Max_Iter)] )
            break % Terminate execution of WHILE .. exit the time evolution loop.

        end    % if(Condition < toler

    end

    if (mod(Cur_Iter,Output_Every)==0) ;  % Output from loop every ...
        %if (Cur_Iter>60) ;  % Output from loop every ...

        rho=sum(f,3); % density
        figure(10); imshow(rho,[0.1 0.9]); title(' rho'); % visualize density evolution
        figure(11); imshow(ux,[ ]); title(' ux' ); % visualize fluid velocity horizontal
        figure(12); imshow(-uy,[ ]); title(' uy' ); % visualize fluid velocity down
        figure(14), imshow(-uyout,[]), title('uyout'); % vis vel flow out
        up=2; % linear section to visualize up from the lower row
        figure(15), hold off, feather(ux(Nr-up,:),uy(Nr-up,:)),
        figure(15), hold on , plot(uy_analy_profile,'r-')
        title('Analytical vs LB calculated, fluid velocity parabolic profile')
        pause(3); % time given to visualize properly

    end % every


   % pause(1);

    
end %  End main time Evolution Loop

% Output & Draw after the end of the time evolution

figure, plot(Cond_path(2:end)); title('convergence path')
%figure, plot(density_path(2:end)); title('density convergence path')
figure, plot( [uy(Nr-up,:)-uy_analy_profile] ); title('difference : LB - Analytical solution')

toc

% Permeability K

K_Darcy_Porous_Sys= (av_vel_int*porosity)/dPdL*Lky_visco ,

K_Analy_2D_Channel=(Width^2)/12


