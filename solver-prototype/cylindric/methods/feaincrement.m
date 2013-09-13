% Procedure for FEA of incompressible axissymetric problems
global g_compressible;
global g_resultspath;
global g_logfile;
global g_jacobi_refinement;

iterations = 100;

if g_compressible 
	compressible = 'compressible';
else 
	compressible = 'incompressible';
end
if g_jacobi_refinement
	jacobi = 'with Jacobi refinement ...';
else
	jacobi = '...';
end
string = sprintf('Starting nonlinear %s job %s',compressible,jacobi);
logger(string);

logger('Constructing elements, geometry and all needed variables');
% Constructing geometry and elements
[nodes,elements,symmetry,boundary] = geometry;
%[nodes,elements,symmetry,boundary] = geometry_test;

% Calculat increment for iterative procedure
dp = max(abs(boundary(:,3)));
if max(boundary(:,3)) == 0
	dp = -dp;
end


[elements_count,corners_count,msize,IDX] = prepare_elements(nodes,elements);


%%
%% Trying to obtain 1st linear solution
%%
logger('Obtaining first iteration from linear job...');
% 1st Linear solution
logger('1st Linear solution - constructing global stiffness matrix...');
[M,P] = construct_global_matrix_lin(msize,nodes,elements,IDX);
logger('1st Linear solution - applying boundary conditions...');
[M,P] = boundary_conditions(M,P,boundary,symmetry,0);
%
% Solve global system
logger('1st Linear solution - solving global system...');
tic;
X = inv(M)*P;
t = toc;
string = sprintf('1st Linear solution - solution time = %f seconds',t);
logger(string);
fname = 'solution_linear1.msh';
string = sprintf('1st Linear solution - exporting to %s file...',fname);
logger(string);
filename = sprintf('%s/%s',g_resultspath,fname);
export_linear(filename,elements,nodes,IDX,X);

logger('Now we have 1st linear solution');

%%
%% Trying to obtain 2nd linear solution
%%
% 2nd Linear solution
logger('2nd Linear solution - constructing global stiffness matrix...');
[M,P] = construct_global_matrix_lin(msize,nodes,elements,IDX);
logger('2nd Linear solution - applying boundary conditions...');
[M,P] = boundary_conditions(M,P,boundary,symmetry,0);
%
% Solve global system
logger('2nd Linear solution - solving global system...');
tic;
X1 = inv(M)*P;
t = toc;
string = sprintf('2nd Linear solution - solution time = %f seconds',t);
logger(string);
fname = 'solution_linear2.msh';
string = sprintf('2nd Linear solution - exporting to %s file...',fname);
logger(string);
filename = sprintf('%s/%s',g_resultspath,fname);
export_linear(filename,elements,nodes,IDX,X1);

logger('Now we have 2nd linear solution');

%%
%% Nonlinear part
%% 

logger('Starting nonlinar part of solution');

% Construct global stiffness matrix
I0 = ones(4,1);
I0(4) = 0;
X0 = X;
X = X1;
% now we have 2 linear solutions: X0, X1

delta = 0.1;

F = zeros(msize,1);

dp = dp/5.0;

[M0,P0] = construct_global_matrix(msize,nodes,elements,IDX,X0);
[M0,P0] = boundary_conditions(M0,P0,boundary,symmetry,0);

for IT = 1:iterations
	logger('\tConstructing global stiffness matrix...');
	tic;
	[M1,P1] = construct_global_matrix(msize,nodes,elements,IDX,X);
	t = toc;
	string = sprintf('\tGlobal matrix construction time = %f seconds',t);
	logger(string);
	[M1,P1] = boundary_conditions(M1,P1,boundary,symmetry,dp*IT);

	f0 = M0*X1-P0;
	f1 = M1*X-P1;
	F = (f0 - f1);% /dp;
	logger('\tConstructing global Jacobi matrix...');
	tic;
	J = construct_global_jacobi(msize,nodes,elements,IDX,X);
	t = toc;
	string = sprintf('\tGlobal Jacobi matrix construction time = %f seconds',t);
	logger(string);
	boundary_jacobi = boundary;
	bnd_siz = size(boundary_jacobi);
	for i = 1:bnd_siz(1)
		if ~boundary_jacobi(i,3) == 0
			boundary_jacobi(i,3) = dp;
		end
	end
	[J,F] = boundary_conditions(J,F,boundary_jacobi,symmetry,0); 
	[dX,flag] = solution(J,F,'\t');	
	
	if not(flag) % solution not converged
		break;
	end
		
	X1 = X;
	X = X + dX;
	M0 = M1;
	P0 = P1;

%	Jacobi solution refinement
    if g_jacobi_refinement
        logger('\t\tJacobi refinement');
		j = 0;
		err = 1;
		X0 = X;
		for i = 1:bnd_siz(1)
			if ~boundary_jacobi(i,3) == 0
				boundary_jacobi(i,3) = 0;
			end
		end

		while err > delta & j < 25 
			logger('\t\tConstructing global Jacobi matrix...');
			AA = zeros(msize,1);
			J = construct_global_jacobi(msize,nodes,elements,IDX,X);
			F = M1*X - P1;
			[J,AA] = boundary_conditions(J,AA,boundary_jacobi,symmetry,0);
			clear AA;
			[dX,flag] = solution(J,-F,'\t\t');	
	    	if not(flag) % solution not converged
				break;
			end

	    	X0 = X;
			X = X + dX;
			err = norm(X-X0);
			string = sprintf('\t\tNorm error = %e, Newton iteration = %d',err,j);
			logger(string);
			% export temporary solution
			fname = sprintf('solution_nonlinear%d_jacobi%d.msh',IT,j);
			filename = sprintf('%s/%s',g_resultspath,fname);
			export_nonlinear(filename,elements,nodes,IDX,X);
			j = j + 1;
		end
	end % Jacobi refinement


	% post-processing
	fname = sprintf('solution_nonlinear%d.msh',IT);
	filename = sprintf('%s/%s',g_resultspath,fname);
	export_nonlinear(filename,elements,nodes,IDX,X);

	workspace_i = sprintf('save %s/iteration%d',g_resultspath,IT);
	eval(workspace_i);

	% error handling	
	err = norm(X - X1);
	if err > 100
		string = sprintf('\tError: norm(X1-X) = %f, calculation stopped',err);
		logger(string);
		break;
	end
	string = sprintf('\tIteration %d done...',IT);
	logger(string);

end
logger('Done!');
