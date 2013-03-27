%................................................................

function outputDisplacementsReactions...
    (displacements,stiffness,GDof,prescribedDof,force)

% output of displacements and reactions in
% tabular form

% GDof: total number of degrees of freedom of 
% the problem

disp('Displacements')
fprintf('node     displacements\n')
% displacements
for jj=1:GDof
    fprintf('%2.0f:        %10.4e\n', jj, displacements(jj))
end

% reactions
F=stiffness*displacements;
reactions=F(prescribedDof)-force(prescribedDof);
disp('Reactions')
fprintf('node     reactions\n')
for jj=1:size(prescribedDof,1)
    fprintf('%2.0f:        %10.4e\n', prescribedDof(jj), reactions(jj))
end
