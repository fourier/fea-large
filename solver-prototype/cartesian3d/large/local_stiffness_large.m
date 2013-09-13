function K = local_stiffness_large(el0,el)
% function K = local_stiffness_large(el0,el)
% calculate local stiffness matrix of element el

    nodes_count = size(el,1);
    dof_count = 3;
    local_size = nodes_count*dof_count;
    K = zeros(local_size,local_size);
    K = integrate5nodes(el,@kmatrix,el0);


function K = kmatrix(el,r,s,t,el0)
    q = create_q_from_elements(el0,el);
    F = graddef(el0,r,s,t,q);
    A = elasticity_matrix_large(F);
    [B,J] = b_matrix(el,r,s,t);
    K = abs(det(J))*B'*A*B;
