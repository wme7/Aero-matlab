Z = peaks;
surf(Z)
axis tight
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;
for k = 1:20
  surf(cos(2*pi*k/20)*Z,Z)
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf) %g443800