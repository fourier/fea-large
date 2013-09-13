function K = construct_k(el,q,r,z)
% Function for constructing k matrix(under integral) 
% by given element el and point (r,z)

B = b_matrix(el,r,z);
B1 = b1_matrix(el,q,r,z);
B2 = b2_matrix(el,q,r,z);

A = elasticity_matrix();

K = (B'+B1'+B2')*A*(B+B1);

