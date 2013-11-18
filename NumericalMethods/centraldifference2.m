
% switch on frame saving
saveframes=1;

% solve the advection equation: du/dt = c*du/dt
% numerical method: RK in time, right differencing
% du_m/dt = c*(u_{m+1}-u_m)/dt

% number of data points
 M = 100;

% periodic interval dimensions (i.e. x \in [a,b))
a = 0; b =  2;

% data point spacing
dx = (b-a)/M;

% coordinates of data points ("nodes")
x = linspace(a, b-dx, M);

% start time
StartTime = 0;
EndTime = 8;

% dt chosen small to make space errors dominant
dt = 0.05*dx;

% choose order of RK integrator
s = 4;

% compute number of time steps
Nsteps = ceil((EndTime-StartTime)/dt);

% modify dt to end exactly at EndTime
dt = (EndTime-StartTime)/Nsteps;

% initial condition 
u = exactsolution(x);

t = StartTime; % set time variable
fram = 1;
for n=1:Nsteps % begin time stepping
    
    utilde = u;  % initiate RK loop

    for rk=s:-1:1 % begin RK sub stages
       
        up1  = [utilde(2:M),utilde(1)];     % utilde_{m+1}
        um1 = [utilde(M),utilde(1:M-1)]; % utilde_{m-1}
        delta0 = 0.5*(up1-um1)/dx;
        
        % update RK variable
        utilde = u+(dt/rk)*delta0;
        
    end % end RK sub stages
  
    u = utilde; % finish RK step
    t  = t+dt;   % update time
    
    if(mod(n,80)==0 |n==Nsteps) % selective plotting
    
        plot(x, u, 'r-*', 'LineWidth', 1); % plot numerical solution
        hold on;

        xplot = linspace(a, b, 1000); % plot exact solution at lots of points
        xmod = mod(xplot+t,b); % coordinate for exact pulse
        uexact = exactsolution(xmod);
        plot(xplot, uexact, 'k-', 'LineWidth', 1);hold off;
        axis([0 2 -0.5 1.5]);
        legend(sprintf('Numerical solution (M=%d,t=%g)', M,t), 'Exact solution');
        drawnow; pause(0.02);
        if(saveframes == 1)        % save frame to file
            if(fram<10)            fname = sprintf('anim_000%d.ppm', fram);
            elseif(fram<100)   fname = sprintf('anim_00%d.ppm', fram);
            elseif(fram<1000) fname = sprintf('anim_0%d.ppm', fram);
            end
            print('-dppm', fname);
            fram = fram+1;
        end
    end
end

% to make the animated gif: I used ImageMagick's convert
% under cygwin:   convert *ppm deltaplus.gif
xmod = mod(x+t,b);
uexact    = exactsolution(xmod);
finalerror = abs(uexact-u); % relative error
