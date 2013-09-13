function R = tvector(el,r,s,t,el0)
% derivatives of shape functions

    q = create_q_from_elements(el0,el);
    % deformation gradient 
    F = graddef(el0,r,s,t,q);
    J = det(F);
    %if J <= 0.0
    %    J
    %    stop_calculation_because_of_gradient
    %end

    [B,jacobi] = b_matrix(el,r,s,t);
    % Cauchy stress tensor
    S = sigma(F);
    % Vector of stresses
    s1 = matrix_to_vector(S);
    R = abs(det(jacobi))*B'*s1;

