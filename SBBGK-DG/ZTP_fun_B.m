function [Zo,To,Po]=ZTP_fun_B(nx,pp,ETo,Ro,Uo)    
%Compute Z, T, P for Boson Gases: IT=-1
Po=zeros(nx,pp);
To=zeros(nx,pp);
Zo=zeros(nx,pp);
for i=1:nx
            for m=1:pp
                ZA = 0.0001;
                ZB = 0.99;
                while (abs(ZA-ZB) > 0.00001)
                    GA12 = 0;
                    GB12 = 0;
                    GA32 = 0;
                    GB32 = 0;
                    for L = 1:50
                            GA12 = GA12 + (ZA^L)/(L^0.5);
                            GB12 = GB12 + (ZB^L)/(L^0.5);
                            GA32 = GA32 + (ZA^L)/(L^1.5);
                            GB32 = GB32 + (ZB^L)/(L^1.5);
                    end
                    PSIA = 2*ETo(i,m) - GA32*(Ro(i,m)/GA12)^3/(2*pi) - Ro(i,m)*Uo(i,m)^2;
                    PSIB = 2*ETo(i,m) - GB32*(Ro(i,m)/GB12)^3/(2*pi) - Ro(i,m)*Uo(i,m)^2;
                    ZC = (ZA + ZB)/2;
                    GC12 = 0;
                    GC32 = 0;
                    GC52 = 0;
                    for L = 1:50
                            GC12 = GC12 + (ZC^L)/(L^0.5);
                            GC32 = GC32 + (ZC^L)/(L^1.5);
                            GC52 = GC52 + (ZC^L)/(L^2.5);
                    end
                    PSIC = 2*ETo(i,m) - GC32*(Ro(i,m)/GC12)^3/(2*pi) - Ro(i,m)*Uo(i,m)^2;
                    
                    if ((PSIA*PSIC) < 0)
                        ZB = ZC;
                    else
                        ZA = ZC;
                    end
                end
                Zo(i,m) = ZC;
                To(i,m) = Ro(i,m)^2 / (pi*GC12^2 );
                Po(i,m) = ETo(i,m) - 0.5 * Ro(i,m) * Uo(i,m)^2;
            end
end
return