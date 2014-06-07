clear all;
clc;
clf;
% coordinates and connectivities 
 x=[0 2 2 0];
 y=[0 0.5 1 1];
 
a = 3;
%mesh number in one side
n=2^(a-1);%(nxn)
% numberElements: number of elements
numberElements=n^2; 
% numberNodes: number of nodes
numberNodes=(n+1)^2;
x1=linspace(x(1),x(2),n+1);
x2=linspace(x(4),x(3),n+1);
y1=linspace(y(1),y(4),n+1);
y2=linspace(y(2),y(3),n+1);
tempY=zeros(n+1,n+1);
tempX=zeros(n+1,n+1);
for i=1:n+1
    tempY(i,:)=linspace(y1(i),y2(i),n+1);
    tempX(:,i)=linspace(x1(i),x2(i),n+1);
end
ct=1;
nodeCoordinates=zeros(numberNodes,2);
for i=1:n+1
    for j=1:n+1
        nodeCoordinates(ct,:)=[tempX(i,j) tempY(i,j)];
        ct=ct+1;
    end
end

elementNodes=zeros(numberElements,4);
for i=1:n
    for j=1:n
        elementNodes(j+n*(i-1),:)=...
            [j+(n+1)*(i-1)...
            j+1+(n+1)*(i-1)...
            j+n+2+(n+1)*(i-1)...
            j+n+1+(n+1)*(i-1)];
    end
end
drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
