function s = stress(nodes,elements,IDX,index,X)
% function s = stress(nodes,elements,IDX,index,X)
q = create_q_from_u(IDX,X,index);
el = element(nodes,elements,index);

a = 0.58541020;
b = 0.13819660;

r = a;
s = b;
t = b;

e = epsilon_lin(el,q,r,s,t);
A = elasticity_matrix();
s = A*e;