function M = elasticity_matrix_large(F)
 M = elasticity_A5(F);
% M = elasticity_neohookean(F);
% M = elasticity_A5_large(F);            
end

function M = elasticity_linear()
% function M = elasticity_matrix()
% Returns Elasticity matrix 6x6 for
% the model A5
% A5: lambda*I1(Epsilon)+2mu*Epsilon

% turn on for elastic model
%if not(true)
    E = 1e9;
    nu = 0.3;
    mul = E/((1+nu)*(1-2*nu));
    M = mul*[1-nu, nu, nu, 0, 0, 0;
             nu, 1-nu, nu, 0, 0, 0;
             nu, nu, 1-nu, 0, 0, 0;
             0, 0,  0,  (1-2*nu)/2., 0, 0;
             0, 0,  0,  0, (1-2*nu)/2., 0;
             0, 0,  0,  0, 0, (1-2*nu)/2.];

end


function M = elasticity_A5(F)
% model A5
    l = 100;
    m = 100;
    M = [l+2*m, l, l, 0, 0, 0;
         l, l+2*m, l, 0, 0, 0;
         l, l, l+2*m, 0, 0, 0;
         0, 0,  0,  m, 0, 0;
         0, 0,  0,  0, m, 0;
         0, 0,  0,  0, 0, m];
    M = M/det(F);
end


function cijkl = cijkl_elasticity(F,i,j,k,l)
% J*(nu*ExE+2*l*D4)
    l = 100;
    m = 100;
    cijkl = det(F)*(l*kroneker_delta(i,j)*kroneker_delta(k,l)+ ...
                    2*m*(kroneker_delta(i,k)*kroneker_delta(j,l)));
end

function M = elasticity_A5_large(F)
    % NOTE: Current implementation is dirty hack!!!
% J*(nu*ExE+2*l*D4)
    M = elasticity_A5();
    C = create_empty_C();
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    %I = tensor_matrix_mapping(i,j);
                    %J = tensor_matrix_mapping(k,l);
                    cijkl = cijkl_elasticity(F,i,j,k,l);
                    %C(i,j,k,l) = elasticity_A5_cijkl(M(I,J),i,j,k,l,F); 
                    C(i,j,k,l) = push_forward_cijkl(cijkl,i,j,k,l,F); 
                end
            end
        end
    end
    for II = 1:6
        for JJ = 1:6
            [i,j] = matrix_tensor_mapping(II);
            [k,l] = matrix_tensor_mapping(JJ);
            M(II,JJ) = C(i,j,k,l);
        end
    end
    M = M/det(F);
end

function cijkl = push_forward_cijkl(old_cijkl,i,j,k,l,F)
    cijkl = old_cijkl;
    for I = 1:3
        for J = 1:3
            for K = 1:3
                for L = 1:3
                    cijkl = cijkl + cijkl*F(i,I)*F(j,J)*F(k,K)*F(l,L);
                end
            end
        end
    end
end


function M = elasticity_neohookean(F)
    M = zeros(6,6);
    C = create_empty_C();
    l = 100;
    m = 100;
    J = det(F);
    l1 = l/J;
    m1 = (m-log(J))/J;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    C(i,j,k,l) = cijkl_neohookean(l1,m1,i,j,k,l);
                end
            end
        end
    end
    for II = 1:6
        for JJ = 1:6
            [i,j] = matrix_tensor_mapping(II);
            [k,l] = matrix_tensor_mapping(JJ);
            M(II,JJ) = C(i,j,k,l);
        end
    end
end

function cijkl = cijkl_neohookean(l1,m1,i,j,k,l)
    cijkl = l1*kroneker_delta(i,j)*kroneker_delta(k,l) + ...
            m1*(kroneker_delta(i,k)*kroneker_delta(j,l) + ...
                kroneker_delta(i,l)*kroneker_delta(j,k));
end

function E4 = create_empty_C()
    E3 = cat(3,zeros(3),zeros(3),zeros(3));
    E4 = cat(4,E3,E3,E3);
end



function I = tensor_matrix_mapping(i,j)
    if i == 1 && j == 1
        I = 1;
    elseif i == 2 && j == 2
        I = 2;
    elseif i == 3 && j == 3
        I = 3;
    elseif (i == 1 && j == 2) || (i == 2 && j == 1)
        I = 4;
    elseif (i == 2 && j == 3) || (i == 3 && j == 2)
        I = 5;
    elseif (i == 1 && j == 3) || (i == 3 && j == 1)
        I = 6;
    end
end

function [i,j] = matrix_tensor_mapping(I)
    if I == 1
        i = 1;
        j = 1;
    elseif I == 2
        i = 2;
        j = 2;
    elseif I == 3
        i = 3;
        j = 3;
    elseif I == 4
        i = 1;
        j = 2;
    elseif I == 5
        i = 2;
        j = 3;
    elseif I == 6
        i = 1;
        j = 3;
    end
end