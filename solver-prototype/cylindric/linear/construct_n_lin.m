function N = construct_n_lin(el,r,z)
% Function for constructing N matrix
% by given element el and point (r,z)

B = b_matrix(el,r,z);
Fp = fp_matrix(el,r,z);
I0 = [1;1;1;0];

N = B'*I0*Fp';

