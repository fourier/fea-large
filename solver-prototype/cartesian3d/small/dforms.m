function [B,J] = dforms(el,r,s,t)
% b = dforms(el,r,s,t)


nodes_count = size(el,1);
dof_count = 3;

M = zeros(dof_count,nodes_count);       % 3x10

% see Zienkiewitz v1, 6th edition, p.146-147
for shape = 1:nodes_count
    for dof = 1:dof_count
        M(dof,shape) = disoform(shape,dof,r,s,t);
    end
end

J = M*el;
% J is a Jacobi matrix of transformation btw local and gloval
% coordinate systems
if abs(det(J)) <= eps
    error('Jacobian of transformation matrix = %f',det(J));
end
J1 = inv3x3(J);

% size_b = nodes_count * dof_count;

B = zeros(dof_count,nodes_count);

%
%   [ dN/dx ]           [ dN/dr ]
%   [ dN/dy ]  = J^-1 * [ dN/ds ]
%   [ dN/dz ]           [ dN/dt ]
%
for i = 1:nodes_count
    p = [disoform(i,1,r,s,t);disoform(i,2,r,s,t);disoform(i,3,r,s,t)];
    b = J1*p;
    B(:,i) = b;
end

