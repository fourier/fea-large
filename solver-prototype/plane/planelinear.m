% Procedure for linear FEA

addpath ('../misc');
% Starting linear task

[nodes,elements,symmetry,boundary] = matlab_geometry6;
%[nodes,elements,symmetry,boundary] = matlab_geometry1;
[elements_count,msize,IDX] = prepare_elements(nodes,elements);

global indd;
indd = 0;

% Construct global stiffness matrix
[M,P] = construct_global_matrix(msize,nodes,elements,IDX);
dumpmatrix('Dump/global.txt',M);
[M,P] = boundary_conditions(M,P,boundary,symmetry);
dumpmatrix('Dump/global_applied.txt',M);
dumpvector('Dump/global_vector.txt',P);


%
% Solve global system
X = inv(M)*P;
dumpvector('solution.txt',X);
% prepare to export to 'msh' file
fname = 'solution_linear.msh';
export_linear(fname,elements,nodes,IDX,X);
