%..............................................................
function outputDisplacementsReactions...
(displacements,stiffness,GDof,prescribedDof)
% output of displacements and reactions in
% tabular form
% GDof: total number of degrees of freedom of
% the problem
% displacements
disp('Displacements')
%displacements=displacements1;
jj=1:GDof; format
[jj' displacements]
% reactions
F=stiffness*displacements;
reactions=F(prescribedDof);
disp('reactions')
[prescribedDof reactions]