for IT1 = 21:200
	nodes = update_with_prescribed(nodes,boundary);

	% Set iterations counter and tolerance
	tolerance = 10000; % set big number to make at least one iteration
	IT2 = 0;

	%
	% Loop over Newton iterations
	% 
	
	while IT2 < newton_max && abs(tolerance) > desired_tolerance
		% Store previous solution
		nodes1 = nodes;

		% Constuct global matrix
		fprintf('Constructing global stiffness matrix and vector...\n');
		M = stiffness(msize,nodes0,nodes,elements,IDX);
		R = residual_force(msize,nodes0,nodes,elements,IDX);  
		[M,R] = apply_bc(M,R,boundary_fixed);
		% Solve global system
		fprintf('Solving global system...\n');
		[X,flag] = solution_qmr(M,-R);
		% Check for convergence
		tolerance = norm(-R - M*X)
		nodes = update_with_solution(nodes,X);
		rtu0 = X'*R;
		tolerance = rtu0
		IT2 = IT2+1;
	end

	fprintf('Exporting data...\n');
	U = create_u_from_nodes(msize,nodes0,nodes,elements,IDX);
%	export_meshviewer('sol.txt',nodes,elements,IDX,X);
	undef_name = sprintf('undeformed%d.msh',IT1);
	def_name = sprintf('deformed%d.msh',IT1);
	export_gmsh_large(undef_name,def_name,nodes0,elements,IDX,U);
end
