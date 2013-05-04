%................................................................

function outputDisplacements...
    (displacements,numberNodes,GDof)

% output of displacements in tabular form

disp('Displacements')
jj=1:numberNodes; format
f=[jj; displacements(1:2:GDof)'; displacements(2:2:GDof)'];
fprintf('Node           UX           UY\n')
fprintf('%4d   %10.4e    %10.4e\n',f)