function Bp = construct_bp(el,q,r,z)
% Function for constructing N matrix
% by given element el and point (r,z)

B = b_matrix(el,r,z);
B1 = b1_matrix(el,q,r,z);
Fp = fp_matrix(el,r,z);

I0 = [1;1;1;0];

Bp = Fp*I0'*(B+B1);

