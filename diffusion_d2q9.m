function diffusion_d2q9
den=0.5; 
%���������Ŀ
nx=128; 
ny=1;
x1=17;
y1=1;
c=1/9;
A=1;
Ts=1;      
step=10;   
while(Ts<=step)
for k=1:9
    for x=1:nx
        for y=1:ny
            iniden=den+A*sin(2*pi*x/nx);
            F(x,y,k)=iniden/9;
            DENSITY(x,y)=iniden;
        end
    end
end
FEQ=F;
Dthe=zeros([1 10]);
Dsim=zeros([1 10]);
D=zeros([1 10]);
nu=length(D);
    omega=1.0/Ts;
    ts=0;
    while (ts<1000) 
    % Propagate
    F(:,:,4)=F([2:nx 1],[ny 1:ny-1],4);
    F(:,:,3)=F(:,[ny 1:ny-1],3);
    F(:,:,2)=F([nx 1:nx-1],[ny 1:ny-1],2);
    F(:,:,5)=F([2:nx 1],:,5);
    F(:,:,1)=F([nx 1:nx-1],:,1);
    F(:,:,6)=F([2:nx 1],[2:ny 1],6);
    F(:,:,7)=F(:,[2:ny 1],7); 
    F(:,:,8)=F([nx 1:nx-1],[2:ny 1],8);
    for ic=1:9;
            FEQ(:,:,ic)=c*DENSITY;
    end
    F=omega*FEQ+(1-omega)*F;
    DENSITY=sum(F,3);
    ts=ts+1;
    end
    %dendiff=DENSITY(x1,y1)-0.5;ln1=log(dendiff);ln2=ln1/(-(2*pi*x1/nx)^2);ln2=log(ln2); Dsim=exp(ln2/ts);
    dendiff=(DENSITY(x1,y1)-0.5)/cos(2*pi*x1/nx)/(40/3/pi);ln1=log(abs(dendiff));Dsim=ln1/(-(2*pi/nx)^2*ts) %��ɢϵ�����ֵ
    %dendiff=(Bave-0.5)/sin(2*pi*x1/nx);ln1=log(abs(dendiff));Dsim=ln1/(-(2*pi/nx)^2*ts)
    Dthe=(Ts-0.5)*2/3   %��ɢϵ������ֵ
    %Dsim=Dsim*Dthe;
    D=abs(Dsim-Dthe)/Dthe      %����ֵ�����ֵ��������
    Dthes(1,Ts)=Dthe;
    Dsims(1,Ts)=Dsim; %#ok<AGROW>
    Ds(1,Ts)=D;
    Ts=Ts+1;
end
%Ts=1:step;
Ts=3:8;
figure(1);
plot(Ts,Dsims(3:8),'o',Ts,Dthes(3:8),'-');grid on
figure(2);
plot(Ts,Ds(3:8),'o');grid on
figure;
[X,Y]=meshgrid(1:nx,1:ny);mesh(X,Y,DENSITY(1:nx,1:ny))