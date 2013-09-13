function R = residual_force(msize,nodes0,nodes,elements,IDX)  
% R = residual_force(msize,nodes,elements,IDX)  
% function for constructing global nodal residual force vector R
% R = T - F
% where T - nodal internal forces, F - nodal external forces
% when no external forces, R = T
%fprintf(['Constructing global vector of residual forces...\n']);
    R = zeros(msize,1);
    elements_count = size(elements,1);

    nodes_count = size(elements,2);
    dof_count = 3;
    local_size = nodes_count*dof_count;


    for i = 1:elements_count
        el = element(nodes,elements,i);
        el0 = element(nodes0,elements,i);
        T = residual(el0,el);
        % distribute in global vector
        for k = 1:local_size
            R(IDX(i,k))=R(IDX(i,k))+T(k);			
        end	
    end


function T = residual(el0,el)
% function T = construct_t(el0,el)
% Function for constuction of internal force vector T
% el - element in current configuration
% q - displacements from initial configuration(for deformation gradient)
% T = int(B^T * Sigma)
    
    T = integrate5nodes(el,@tvector,el0);
    
function r = tvector(el,r,s,t,el0)
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
    s = matrix_to_vector(S);
    r = abs(det(jacobi))*B'*s;
