% Roe average Euler solver
% by Yun-da Tsai 2012/12/21  Sanyi, Mioli

% Roe method ref:A course in Computational Fluid Dynamics


% Parameters Definition
% system eq: setting q
ncells=80;
rho=zeros(1,80);
rho_u=zeros(1,80);
rho_E=zeros(1,80);
tf = 1.7; %Finial time
cfl = 0.08; % CFL number
% Space Domain
xmin =0.0;
xmax =1.0;
dx = (xmax-xmin)/ncells;
gamma=1.4;
Cv=1.4182;%Cp
R=8.314;


% Initial condition 

for i=1:ncells
    if(i<=(ncells/2))
        q_Initial(1,i)=1;  %Density
        q_Initial(2,i)=0;  %Velocity
        q_Initial(3,i)=1;  %Pressure
    else
        q_Initial(1,i)=0.125;
        q_Initial(2,i)=0;
        q_Initial(3,i)=0.1;
    end
end

% BC
q_Initial_previous(1:3,1) = q_Initial(1:3,1);
q_Initial_later(1:3,1) = q_Initial(1:3,80);

% Use initial condition to representthe matrix q
% q(rho/rho*u/rho*e)
q(1,:)=q_Initial(1,:);
for i=1:ncells
 q(2,i) =q_Initial(2,i)*q_Initial(1,i);    
 q(3,i)=q_Initial(3,i)/((gamma-1)+0.5*q_Initial(1,i)*(q_Initial(2,i))^2);
end

% Calculate the Roe Average
% D_bar 
for i=1:79
      D_bar(1,i)=sqrt(q_Initial(1,i+1)/q_Initial(1,i));    
end
% U_bar
for j=1:79
    Temp_1=q_Initial(2,i)+D_bar(1,i)*q_Initial(2,i+1);
    u_bar(1,i)=Temp_1/(1+D_bar(1,i));
end

% rho_bar

for i=1:79    

    rho_bar(1,i)=  q_Initial(1,i)*((0.5*(1+D_bar(1,i))).^2);
end

% h_bar
for i=1:79
   h_bar(1,i)= (((Cv*q_Initial(3,i))/ (q_Initial(1,i)*R))+D_bar(1,i)*Cv*(q_Initial(3,i+1)/(q_Initial(1,i+1)*R)))/(1+D_bar(1,i));
end
% a_bar^2
for i=1:79    
    a_bar_sqaure(1,i)=(gamma-1)*(h_bar(1,i)-0.5*(u_bar(1,i).^2));
end
% p_bar
for i=1:79
    p_bar(1,i)=(1/gamma)*rho_bar(1,i)*a_bar_sqaure(1,i);
end
% e_bar & E_bar
for  i=1:79
E_bar(1,i)=rho_bar(1,i)*h_bar(1,i)-p_bar(1,i);    
e_bar(1,i)=a_bar_sqaure(1,i)/(gamma*(gamma-1));
end

% Coefficient Matrix  3X3 A(u)

for i=1:79
              Temp_10=a_bar_sqaure(1,i)/(gamma-1);
        A=[0                                                                             1                                                                                                                0;
              ((gamma-3)*(u_bar(1,i)).^2)/2                                      -(gamma-3)*u_bar(1,i)                                                           gamma-1
              u_bar(1,i)*(((gamma-2)*(u_bar(1,i)).^2)/2-Temp_10)   (3-2*gamma* (u_bar(1,i).^2))/2+Temp_10                   gamma* u_bar(1,i)];          
               [s,t]= eig(A');
               eval(['s' num2str(i) '=s']);
end
% Eigenvalue
eigenvalue_1=u_bar;
eigenvalue_2=u_bar-sqrt(a_bar_sqaure);
eigenvalue_3=u_bar+sqrt(a_bar_sqaure);
eigenvalue=[eigenvalue_1 ; eigenvalue_2;eigenvalue_3];

for i=1:79
    Temp_5=sqrt(a_bar_sqaure);
    Temp_6(1,i)=u_bar(1,i)*Temp_5(1,i);
end

eigenvector_1=[ones(1,79) ; u_bar ; 0.5*(u_bar.^2)];
eigenvector_2=[ones(1,79); eigenvalue_2; h_bar-Temp_6];
eigenvector_3=[ones(1,79); eigenvalue_2 ; h_bar+Temp_6];

% eigenvector=[eigenvector_1 ; eigenvector_2; eigenvector_3];

for i=1:79
      delta(1,i)=q_Initial(1,i+1)-q_Initial(1,i);
      delta(2,i)=q_Initial(2,i+1)-q_Initial(2,i);
      delta(3,i)=q_Initial(3,i+1)-q_Initial(3,i);
end

for i=1:79    
    	Final= eval(['s' num2str(i)]); 	       
        alpha(1,i)= Final(1:3,1)'*delta(1:3,i);
        alpha(2,i)= Final(1:3,2)'*delta(1:3,i);
        alpha(3,i)= Final(1:3,3)'*delta(1:3,i);
end



% for i=1:timeStep
%     if (i~=tf)
%         
%     end    
% end
% 
for i=1:79
     lamda_plus(1,i)=0.5*(eigenvalue(1,i)+abs(eigenvalue(1,i)));
     lamda_plus(2,i)=0.5*(eigenvalue(2,i)+abs(eigenvalue(2,i)));
     lamda_plus(3,i)=0.5*(eigenvalue(3,i)+abs(eigenvalue(3,i)));
     
     lamda_minus(1,i)=0.5*(eigenvalue(1,i)-abs(eigenvalue(1,i)));
     lamda_minus(2,i)=0.5*(eigenvalue(2,i)-abs(eigenvalue(2,i)));
     lamda_minus(3,i)=0.5*(eigenvalue(3,i)-abs(eigenvalue(3,i)));
end


% 
        t = 0;                           % Initialize the current time.
%        nsteps = 0;                        %!Initialize the number of time steps.
%        time_step=0; 
%        for itime = 1:50000 %50000 is large enough to reach tf=1.7.
%         if(t==tf)
%             break;
%         end        
%         third=subroutine();
%         if((t+dt)>tf)
%             dt=tf-t;
%             t=t+dt;
%         end        
%        end
  [R C]=find(eigenvalue==max(max(eigenvalue)));  
   dt = dx*cfl/eigenvalue(R(1),C(1)); 
       
% Compute next time step
for i=1:(tf/dt)
%     u_next = u - dtdx*(F_right - F_left) ...
 %                     + (dt/r_time)*(u_eq-u);              
       
%        for k=1:39
%                   q_next(1:3,k)= q(1:3,k)-(dt/dx)*(lamda_plus(1,k)*alpha(1,k)*eigenvector_1(1:3,k)...
%                       +lamda_plus(2,k)*alpha(2,k)*eigenvector_2(1:3,k)+lamda_plus(3,k)*alpha(3,k)*eigenvector_3(1:3,k));
%        end
       for j=1:79
                  q_next(1:3,j)= q(1:3,j)-(dt/dx)*(lamda_minus(1,j)*alpha(1,j)*eigenvector_1(1:3,j)...
                      +lamda_minus(2,j)*alpha(2,j)*eigenvector_2(1:3,j)+lamda_minus(3,j)*alpha(3,j)*eigenvector_3(1:3,j)...
                      +lamda_plus(1,j)*alpha(1,j)*eigenvector_1(1:3,j)+lamda_plus(2,j)*alpha(2,j)*eigenvector_2(1:3,j)...
                      +lamda_plus(3,j)*alpha(3,j)*eigenvector_3(1:3,j));
       end   
       eval(['qSimulation' num2str(i) '=q_next']);
end

distance=1:79;

plot (distance (1,distance),(qSimulation444 (1,distance)),'k->');
