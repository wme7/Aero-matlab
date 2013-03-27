Fs = 100;
   dt = 1/Fs;
   StartTime = -5;
   StopTime = 15;
   t = StartTime:dt:StopTime-dt;
   x = (t>1) - (t>2) + (t<2) - (t<3) + (t>3) - (t>4);
   figure;
   stairs(t,x);
   ylim([-1.2 1.2]);