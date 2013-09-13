function e = epsilon_lin(el,q,r,s,t)
% function e = epsilon_lin(el,q,r,s,t)

B = b_matrix(el,r,s,t);
e = B*q;

