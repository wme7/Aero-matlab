function z=adiscontnuous_function(x,y)
z=((abs(y-0.6)+abs(x))<0.4)+(abs(x)<0.05).*(abs(y-0.3)<0.6);