function M = stiffness(msize,nodes0,nodes,elements,IDX)
% M = stiffness(msize,nodes0,nodes,elements,IDX)
% Function for constructing global nonlinear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix
fprintf(['Constructing global stiffness matrix...\n']);

M = zeros(msize);
M = sparse(M);

elements_count = size(elements,1);

for i = 1:elements_count
    el = element(nodes,elements,i);
    el0 = element(nodes0,elements,i);

    Kc = local_stiffness_large(el0,el);
    Ks = local_sigma(el0,el);
    
    K = Kc + Ks;

    local_size = size(Kc,1);
    for k = 1:local_size
        for l = 1:local_size
            M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + K(k,l);
        end
    end
end


