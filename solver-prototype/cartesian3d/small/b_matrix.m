function [B,J] = b_matrix(el,r,s,t)
% function B = b_matrix(el,r,s,t)
% Calculate B matrix - matrix of gradients
% by given element el and local coordinates r,s,t

nodes_count = size(el,1);
dof_count = 3;
[b,J] = dforms(el,r,s,t);


B = zeros(6,nodes_count*dof_count);
for i = 1:nodes_count
  col1 = i*3-2;
  col2 = i*3-1;
  col3 = i*3;
  B(1,col1) = b(1,i);
  B(2,col2) = b(2,i);
  B(3,col3) = b(3,i);
  B(4,col1) = b(2,i);
  B(4,col2) = b(1,i);
  B(5,col2) = b(3,i);
  B(5,col3) = b(2,i);
  B(6,col1) = b(3,i);
  B(6,col3) = b(1,i);
end

