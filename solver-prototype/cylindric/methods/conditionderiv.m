%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure for FEA of incompressible axissymetric problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Import global variables
% 

global g_compressible;
global g_resultspath;
global g_logfile;
global g_linesearch;

% 
% Set calculation parameters
% 

% Global increment iterations
iterations = 50;
% Maximum number of Newton iterations
newton_max = 50;
% Desired tolerance for Newton iterations
desired_tolerance = 1e-5;
% Maximum number of line searches
linesearch_max = 5;

%
% Form log header
%

if g_compressible 
	compressible = 'compressible';
else 
	compressible = 'incompressible';
end
if g_linesearch
	linesearch = 'with line search procedure ...';
else
	linesearch = '...';
end
string = sprintf('Starting nonlinear %s job %s',compressible,linesearch);
logger(string);

%
% Importing and constructing geometry and elements
% 

logger('Constructing elements, geometry and all needed variables');
[nodes,elements,symmetry,boundary] = matlab_geometry;
%[nodes,elements,symmetry,boundary] = geometry;
boundary = convert_boundary(boundary);
[elements_count,corners_count,msize,IDX] = prepare_elements(nodes,elements);

% 
% Solve linear task - just to export it to a file.
% Solution from linear task will not be used in method
% 
logger('Starting linear part of solution');
[M,P] = construct_global_matrix_lin(msize,nodes,elements,IDX);
M0 = M;
P0 = P;
[M,P] = boundary_conditions(M,P,boundary,symmetry,0);
X0 = inv(M)*P;
fname = 'linear.msh';
filename = sprintf('%s/%s',g_resultspath,fname);
export_linear(filename,elements,nodes,IDX,X0);

%%
%% Nonlinear part
%% 

logger('Starting nonlinar part of solution');

% 
% Preserve initial geometry and set solution to 0
%

nodes0 = nodes;
X = zeros(msize,1);

%
% Loop over increments
%

for IT1 = 1:iterations
	
	% Update nodes with prescribed displacements	
	nodes = update_with_prescribed(nodes,boundary);

	% Set iterations counter and tolerance
	tolerance = 10000; % set big number to make at least one iteration
	IT2 = 0;
	rtu0 = 0;

	%
	% Loop over Newton iterations
	% 
	
	while IT2 < newton_max && tolerance > desired_tolerance
		% Store previous solution
		nodes1 = nodes;

		% Constuct global matrix
		logger('\t\tConstructing global stiffness matrix...');
		tic;
		M = construct_global_matrix_cond(msize,nodes0,nodes,elements,IDX);
		P = construct_global_residual_force(msize,nodes0,nodes,elements,IDX);
		M0 = M;
		P0 = P;
		t=toc;
		string=sprintf('\t\tGlobal matrix construction time = %f seconds',t);
		logger(string);
		% Apply boundary conditions in form of displacements
		% with fixed nodes and symmetry
		[M,P] = boundary_conditions_fixed(M,P,boundary,symmetry);
		% Solve a system
		[X,flag] = solution(M,-P,'\t\t');	
		% Check for convergence
		tolerance = norm(-P - M*X);
		if not(flag) % solution method not converged
			if tolerance > desired_tolerance
				string = sprintf('\t\tSLAE method not converged. ||-P-M*X|| = %e',...
					tolerance);
				logger(string); 
				break;
			end	
		end
		string = sprintf('\t\t||-P-M*X|| = %e',tolerance);
		logger(string);
		string = '\t\tUpdate nodes with solution';logger(string);
		nodes = update_with_solution(nodes,elements,IDX,X);

		% get and store a residual value
		rtu0 = X'*P;
		tolerance = abs(rtu0); 
		string=sprintf('\t\tNewton iteration %d finished with tolerance %e', ...
			IT2,tolerance);logger(string);	

		if tolerance <= desired_tolerance % already converged
			break;													% no need to iterate more
		end	
		string=sprintf('\t\tNewton iteration %d finished with tolerance %e', ...
			IT2,tolerance);logger(string);	

		% Not converged at desired tolerance
		
		%
		% Line search method
		%
		if g_linesearch			
			rho = .5;
			eta0 = 0.0 ; eta = 1.d0;
			rtu = rtu0 * rho * 2;
			IT3 = 1;
			X1 = X;
			while IT3 <= linesearch_max && abs(rtu) > abs(rtu0 * rho)
				string = sprintf('\t\t\tLine search iteration %d',...
					IT3);
				logger(string);
				X1 = (eta-eta0)*X;
				nodes = update_with_solution(nodes1,elements,IDX,X1);
				R = construct_global_residual_force(msize,nodes0,nodes,elements,IDX);
				[rtu,eta,eta0] = line_search(R,X1,eta0,eta,rtu0,rtu);
				string = sprintf('\t\t\trtu = %f, rtu0*rho = %f, eta-eta0 = %f ',...
					rtu,rtu0*rho,eta-eta0);
				logger(string);

				IT3 = IT3 + 1;
			end
		end

		% Increment iteration counter
		IT2 = IT2+1;
	end		
	
	%
	% Check for convergence
	%

	if IT2 == newton_max && tolerance > desired_tolerance
		string = sprintf('\tNewton method not converged with tolerance %e', ...
			tolerance); 			
		logger(string);	
		break;			
	else
		string = sprintf('\tNewton method converged with tolerance %e',tolerance);
		string = sprintf('%s on iteration %d',string,IT2);
		logger(string);		
	end	

	%
	% Do some post-processing
	% 

	% Export to msh file
	fname = sprintf('nonlinear%d.msh',IT1);
	filename = sprintf('%s/%s',g_resultspath,fname);
	export_condder(filename,elements,nodes0,nodes,boundary,IDX,msize);
	% Export MATLAB workspace
	workspace_i = sprintf('save %s/iteration%d',g_resultspath,IT1);
	eval(workspace_i);
	% Finalize iteration
	string = sprintf('\tIteration %d done...',IT1);
	logger(string);
end
% well done
logger('Done!');
