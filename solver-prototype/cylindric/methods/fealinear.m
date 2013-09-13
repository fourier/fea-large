% Procedure for linear FEA of incompressible axissymetric problems
global g_compressible;
global g_resultspath;

% Starting linear task
% Drop to log task description
if g_compressible 
	compressible = 'compressible';
else 
	compressible = 'incompressible';
end
string = sprintf('Starting linear %s job',compressible);
logger(string);

logger('Preparing geometry...');

%[nodes,elements,symmetry,boundary] = geometry;
[nodes,elements,symmetry,boundary] = matlab_geometry;

[elements_count,corners_count,msize,IDX] = prepare_elements(nodes,elements);

% Construct global stiffness matrix
logger('Constructing global stiffness matrix...');
[M,P] = construct_global_matrix_lin(msize,nodes,elements,IDX);
dumpmatrix('results/global_matrix.txt',M);
M1 = M;
% apply boundary conditions
logger('Applying boundary conditions...');
[M,P] = boundary_conditions(M,P,boundary,symmetry,0);
dumpmatrix('results/global_matrix_applied.txt',M);
dumpvector('results/global_vector.txt',P);

%
% Solve global system
logger('Solving global system...');
tic;
X = inv(M)*P;
t = toc;
string = sprintf('Solution time = %f seconds',t);
logger(string);
R = M*X - P;
dumpvector('results/residual.txt',R);

% prepare to export to 'msh' file
fname = 'solution_linear.msh';

string = sprintf('Exporting to %s file...',fname);
logger(string);

filename = sprintf('%s/%s',g_resultspath,fname);

export_linear(filename,elements,nodes,IDX,X);

logger('Done!');
