clear
p=3;%polynomial degree
M=6;%number of refinement
n=4;
for l=1:M
 np(l)=l
 %number of intervals
 [err_loc(l,:), jmp(l,:)]=dg_ode2(n,np(l));
end
for l=1:M
    for k=1:n
    ratio(l,k)=err_loc(l,k)/jmp(l,k);
end
end
%err_lvl
%err_lvl
%ord
%T=log10(n);
%[c,d]=polyfit(T(3:M),log10(err_lvl(3:M)),1)
%[cn,dn]=polyfit(T(3:M-1),log10(err_node(3:M-1)),1)
%set(axes,'fontsize',18)
%hold on
%plot(T,log10(err_lvl),'s');
%plot(T,log10(err_node),'o');
%err_fit=c(1)*T+c(2);
%errn_fit=cn(1)*T+cn(2);
%plot(T,err_fit);
%plot(T,errn_fit,'-.');
%xlabel('-log_{10} N');
%ylabel('log_{10} err');
%hold off
%%fid = fopen('order.txt','w');
%%sprintf('N  h   L2  order')
%%sprintf( '%d %6.4f %6.4f  %6.3f  %6.4f %6.4f',n(1),h(1),L1(1),0,L2(1),0)
