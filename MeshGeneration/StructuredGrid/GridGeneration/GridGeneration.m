% To demonstrate the grid generation using analytic coordinate systems
%{
 Author : Siva Srinivas Kolukula                                
          Senior Research Fellow                                
          Structural Mechanics Laboratory                       
          Indira Gandhi Center for Atomic Research              
          India                                                 
 E-mail : allwayzitzme@gmail.com                                         
          http://sites.google.com/site/kolukulasivasrinivas/                 
%}

% Reference: Fundametnals of Grid Generation - Knupp, Steinberg

clc 
clear all
% number of discretizations along xi and eta axis
m = 20 ;
n = 20 ;

% discretize along xi and eta axis
xi = linspace(0.,1,m) ;
eta = linspace(0.,1.,n) ;

Analytic = [{'PolarCoordinates'} ;{'ParabolicCylinderCoordinates'} ;
            {'EllipticCylinderCoordinates'}; {'Horseshoe'} ;
             {'ModifiedHorseshoe'} ;{'BipolarCoordinates'}] ;
         
for grid = 1:length(Analytic) 
    
    % Initialize matrices in x and y axis
    X = zeros(m,n) ;
    Y = zeros(m,n) ;
    % Run a loop alon xi and eta axis to get x,y 
    for i = 1:m
        for j = 1:n
            Xi = xi(i) ;
            Eta = eta(j) ;
            [x y] = feval(Analytic{grid},Xi,Eta) ;

            X(i,j) = x ;
            Y(i,j) = y ;
    
        end
    end
    % To plot grid obtained
    plotgrid(X,Y)
    title(Analytic{grid},'color','b')
    disp('press any key to plot next grid')
    pause
end

