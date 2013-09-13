function h = construct_h(el,q,r,z)
% Function for constructing H vector 
% by given element el,solution q and point (r,z)

Fp = fp_matrix(el,r,z);
e = epsilon_nonlin(el,q,r,z);

h = Fp*h_function(e);

