% animation with drawnow
a = [1:100];
for i=1:100,
 plot([1:i], a(1:i));
 drawnow
end