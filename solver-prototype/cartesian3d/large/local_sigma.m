function K = local_sigma(el0,el)
% K = local_sigma(el0,el)
    
    nodes_count = size(el,1);
    dof_count = 3;
    local_size = nodes_count*dof_count;
    K = zeros(local_size,local_size);
    K = integrate5nodes(el,@ksigma,el0);

function ks = ksigma(el,r,s,t,el0)
    nodes_count = size(el,1);
    dof_count = 3;
    local_size = nodes_count*dof_count;

    ks = zeros(local_size,local_size);
    q = create_q_from_elements(el0,el);
    % deformation gradient
    F = graddef(el0,r,s,t,q);
    J = det(F);
    %if J <= 0.0
    %	J
    %	stop_calculation_because_of_gradient
    %end	
    
    Sigma = sigma(F);
    E = eye(3,3);

    [B,jacobi] = dforms(el,r,s,t);

    for k = 1:nodes_count
        for l = 1:nodes_count
            Bk = B(:,k);
            Bl = B(:,l);
            Hkl = abs(det(jacobi))*Bk'*Sigma*Bl*E;
            ks(k*dof_count-2:k*dof_count,l*dof_count-2:l*dof_count) = Hkl;
        end
    end	
    