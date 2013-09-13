% Procedure for linear FEA

%addpath ('../misc');
% Starting linear task
fprintf('Loading geometry...\n');
[nodes,elements,boundary,boundary_forces] = brick_new3;
boundary = update_boundary(nodes,boundary);

IDX = index_matrix(elements);

% Number of elements
elements_count = size(elements,1);
% Number of nodes per element
nodes_per_element = size(elements,2);
dof = size(nodes,2);
msize = size(nodes,1)*dof;

fprintf('Constructing global stiffness matrix and vector...\n');
% Construct global stiffness matrix
% msize-by-msize all zero sparse matrix
M = sparse(msize,msize);
% right-part vector
P = zeros(msize,1);
for i = 1:elements_count
    el = element(nodes,elements,i);
    K = local_stiffness(el);
    local_size = size(K,1);
    for k = 1:local_size
        for l = 1:local_size
            M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + K(k,l);
        end
    end
end
fprintf(['Condition number of M before applying boundary conditions: ' ...
         '%e\n'],condest(M));

%dumpmatrix('Dump/global.txt',M);
[M,P] = apply_bc2(M,P,boundary);
%dumpmatrix('Dump/global_applied.txt',M);
%dumpvector('Dump/global_vector.txt',P);
fprintf(['Condition number of M after applying boundary conditions: ' ...
         '%e\n'],condest(M));

%
% Solve global system
fprintf('Solving global system...\n');
%X = inv(M)*P;
%return
[L,U] = luinc(M,1e-14);
[X,flag,relres,iter] = pcg(M,P,1e-15,10000,L,U);
check_solution(nodes,elements,boundary,X);

fprintf('Exporting data...\n');
export_meshviewer('solution.txt',nodes,elements,IDX,X);
export_gmsh('undeformed.msh','deformed.msh',nodes,elements,IDX,X);
export_stresses_small('stresses.txt',nodes,elements,IDX,X);
