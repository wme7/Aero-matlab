 function [filter_sigma]=filter_profile(PolyDeg,filter_order, CutOff)
 % declare local arguments
%    integer :: i,j
%    real(kind=8) :: machine_zero, tmp
filter_sigma=zeros(1,1,PolyDeg+1);
    machine_zero=1.0E-10;
    filter_alpha=-log(machine_zero);
    
    DegCutOff=round(CutOff*PolyDeg)
    
     
    for i=0:PolyDeg
       if ( i > DegCutOff ) 
          tmp=((i-DegCutOff)/(PolyDeg-DegCutOff))^filter_order;
          filter_sigma(i+1)=exp(-filter_alpha*tmp);
       else
          filter_sigma(i+1)=1.d0;
       end
       %         filter_sigma(i)=1.d0
      end

%     for i=0:PolyDeg
%        if ( filter_sigma(i+1) <= 1.0E-14 )  
%            filter_sigma(i+1)=0.d0;
%        end
%     end

    return
