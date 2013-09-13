function K = construct_k(el,x,y)

global indd;

B = b_matrix(el,x,y);
A = elasticity_matrix();

fname = sprintf('Dump/B%d.txt',indd);
dumpmatrix(fname,B);

K = B'*A*B;

indd = indd + 1;