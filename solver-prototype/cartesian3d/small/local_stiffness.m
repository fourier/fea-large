function K = local_stiffness(el)
% function K = local_stiffness(el)
% calculate local stiffness matrix of element el
    
    nodes_count = size(el,1);
    dof_count = 3;
    local_size = nodes_count*dof_count;
    K = zeros(local_size,local_size);
    % the choice of number of integration point % described in 
    % Zienkiewitz v1, 6th edition, p.164-168
    K = integrate5nodes(el,@kmatrix);
    %S = volume(el);
    %K = S*K;
end

function k = kmatrix(el,r,s,t)
    A = elasticity_matrix();
    [B,J] = b_matrix(el,r,s,t);
    k = abs(det(J))*B'*A*B;
end
