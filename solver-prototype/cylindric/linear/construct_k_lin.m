function K = construct_k_lin(el,r,z)
% function K = construct_k_lin(el,r,z)
% Function for constructing k matrix(under integral) 
% by given element el and point (r,z)

% In case of small deformations we assume deformation
% gradient as 1 
F = eye(3,3);

fprintf('%.8f %.8f\n',r,z);
B = b_matrix(el,r,z);
A = elasticity_matrix(F);

K = B'*A*B;

