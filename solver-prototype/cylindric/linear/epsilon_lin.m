function e = epsilon_lin(el,q,r,z)
% Function for calculating linear deformations by given
% element el and vector of displacements q, return vector e, where
% e1 => epsilon_{rr}, e2 => epsilon_{zz}, e3 => epsilon_{ff},
% e4 => epsilon_{rz}, 
% by given point r,z.

B = b_matrix(el,r,z);
e = B*q;

