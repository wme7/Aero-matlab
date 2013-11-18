function [p, t]=fixmesh(p,t)

% FIXMESH deletes double points produced in distmesh2d
%   [p, t] = fixmesh(p,t) deletes double points produced by distmesh2d. It
%   also adjusts the triangulation matrix t accordingly. Points in t are
%   also permuted to produce counter clockwise oriented triangles.
%
%   M. Truffer, Oct. 2004
%
%   see also: DISTMESH2D

% check for doubles

Np=size(p,1); dble=zeros(1,Np);
for i=1:Np-1
    for j=i+1:Np 
        if norm(p(i,:)-p(j,:))==0, 
            dble(j)=1;
        end 
    end 
end

% delete the second point of any doubles (t only seems to refer to the
% first point)

i=1;
while i<=size(dble,2)
    if dble(i) 
        p(i,:)=[]; 
        [I J]=find(t>i);
        % all points > i need to be reduced by 1.
        for k=1:size(I,1)
            t(I(k),J(k))=t(I(k),J(k))-1; 
        end
        dble(i)=[]; 
    end
    i=i+1; 
end

% check for the order of nodes in a triangle. If clockwise permute node 2
% and 3

for q=1:size(t,1)
    j=t(q,1); k=t(q,2); l=t(q,3);
    b2=p(l,2)-p(j,2); b3=p(j,2)-p(k,2); 
    c2=-p(l,1)+p(j,1); c3=-p(j,1)+p(k,1);
    twoA=b2*c3-c2*b3;
    if twoA<0
        t(q,2)=l; t(q,3)=k;  
    end 
end