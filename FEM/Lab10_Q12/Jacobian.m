% ............................................................. 

    function [JacobianMatrix,invJacobian,XYDerivatives]=...
        Jacobian(nodeCoordinates,naturalDerivatives)

    % JacobianMatrix    : Jacobian matrix
    % invJacobian : inverse of Jacobian Matrix
    % XYDerivatives  : derivatives w.r.t. x and y
    % naturalDerivatives  : derivatives w.r.t. xi and eta
    % nodeCoordinates  : nodal coordinates at element level
   
    JacobianMatrix=nodeCoordinates'*naturalDerivatives;                   
    invJacobian=inv(JacobianMatrix);
    XYDerivatives=naturalDerivatives*invJacobian;

    end % end function Jacobian
