% Procedure for nonlinear FEA

addpath ('../small');
addpath ('../');
addpath('../../misc');
% Starting task
fprintf('Loading geometry...\n');
[nodes,elements,boundary0] = brick_new2;
% update boundary to set right-side counditions only by Y axis
boundary = update_boundary(nodes,boundary0);

% create an index matrix
IDX = index_matrix(elements);

% Global increment iterations
iterations = 100;
% Desired tolerance for Newton iterations
desired_tolerance = 1e-8;
% Maximum number of line searches
linesearch_max = 0;
% parameter for linesearch method
rho = 0.5;
% Use modified or ordinary Newthon method
modified_newton = true;
% store initial geometry
nodes0 = nodes;


% Number of elements
elements_count = size(elements,1);
% Number of nodes per element
nodes_per_element = size(elements,2);
% Number of degrees of freedom - this number 
% constructed 'formally' from nodes
% indeed it is always 3
dof = size(nodes,2);
% size of global matrix and vector
msize = size(nodes,1)*dof;

% Prepare global stiffness matrix and vector
M = zeros(msize,msize);
M = sparse(M);
R = zeros(msize,1);

% constucting fixed BCs
fixed_boundary = boundary;%boundary_fixed(boundary);
fixed_boundary(:,3) = 0.0;

for IT1 = 1:iterations
    % recalculate nodes to current configuration
    nodes = update_with_prescribed(nodes0,nodes,boundary,IT1);
    % Set inner iterations counter and tolerance
    IT2 = 0;
    M = stiffness(msize,nodes0,nodes,elements,IDX);
    R = - residual_force(msize,nodes0,nodes,elements,IDX);  
    tau = min(diag(M));
    Xprev = zeros(msize,1);
    % Loop over Newton iterations
    tolerance = norm(R);
    while abs(tolerance) > desired_tolerance
        IT2 = IT2+1;
        tic;
        fprintf('Iteration: load %d, internal loop %d\n',IT1,IT2);
        % Store previous solution
        nodes1 = nodes;
        % Constuct global matrix
        % when not using modified Newton method recalc stiffness
        if (~ modified_newton) && IT2 > 1
            M = stiffness(msize,nodes0,nodes,elements,IDX);
        end
        R = - residual_force(msize,nodes0,nodes,elements,IDX);  
        [M,R] = apply_bc2(M,R,fixed_boundary);
        % Solve global system
        fprintf('solving global system...\n');
        X = M\R;
        fprintf('|M*X-R| = %e\n',norm(M*X-R));
        fprintf('inv(M)*X = %e\n',diag(inv(M))*
        % Check for convergence
        nodes = update_with_solution(nodes1,X);
        R1 = residual_force(msize,nodes0,nodes,elements,IDX);
        tolerance = X'*R;
        fprintf('Tolerance(<X,R>) = %e\n',tolerance);
        fprintf('|X-Xprev| = %e\n', norm(X-Xprev));
        Xprev = X;
        eta = 1; % multiplier for the solution X
        % line search
        if linesearch_max > 0
            a = 0.5;
            b = 1;
            tau = (sqrt(5)-1)/2.;
            x1 = a;
            x2 = b;
            for linesearch_iter = 1:linesearch_max
                x1 = b - tau*(b-a);
                x2 = a + tau*(b-a);
                % calculate f(a)
                nodes = update_with_solution(nodes1,X*x1);
                R = -residual_force(msize,nodes0,nodes,elements,IDX);
                f1 = abs(x1*X'*R);
                fprintf('f(x1 = %f)=%e\n',x1,f1);
                % calculate f(b)
                nodes = update_with_solution(nodes1,X*x2);
                r = -residual_force(msize,nodes0,nodes,elements,IDX);
                f2 = abs(x2*X'*R);
                fprintf('f(x2 = %f)=%e\n',x2,f2);
                % test
                if f1 > f2
                    a = x1;
                else
                    b = x2;
                end
                if abs(tolerance) < f1 && abs(tolerance) < f2
                    eta = 1;
                    break
                end
                fprintf('a = %f, b = %f\n',a,b);
                eta = (x1 + x2)/2.;
            end
        end % end of line search algorithm
        nodes = update_with_solution(nodes1,X*eta);

        % update nodes to iteration results, current configuration
        %nodes = update_with_solution(nodes,X);
        it_time = toc;fprintf('Iteration time: %f seconds\n',it_time);
    end
    fprintf('Exporting data...\n');
    U = create_u_from_nodes(msize,nodes0,nodes,elements,IDX);
    undef_name = sprintf('undeformed%d.msh',IT1);
    def_name = sprintf('deformed%d.msh',IT1);
    export_gmsh_large(undef_name,def_name,nodes0,elements,IDX,U);
    stresses_name = sprintf('stresses%d.txt',IT1);
    export_stresses(stresses_name,nodes0,elements,IDX,U);

end
