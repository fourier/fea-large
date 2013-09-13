function m = boundary_conditions_jacobi(M,boundary,symmetry)
% Apply boundary conditions to the Jacobi matrix M and 
% right-side vector P
% Parameters:
% M - global Jacobi matrix
% boundary - matrix of boundary conditions
% symmetry - matrix of nodes under symmentry condition(r- fixed)
% used in iterative methods for solution of nonlinear tasks

m = M;

boundary_count = size(boundary);
boundary_count = boundary_count(1);
msize = length(m);

% apply r-coordinate conditions
for i = 1:boundary_count
	index = boundary(i,1)*2-1; % r-coordinate
	tmp = m(index,index);
	for j = 1:msize
		m(index,j) = 0;
		m(j,index) = 0;
	end
	m(index,index) = tmp;
end
% apply z-coordinate conditions
for i = 1:boundary_count
   	tmp = m(index,index);
	index = boundary(i,1)*2; % z-coordinate
	tmp = m(index,index);
	for j = 1:msize
		m(index,j) = 0;
		m(j,index) = 0;
	end
	m(index,index) = tmp;
end
% apply symmetry conditions
sym_count = size(symmetry);
for i=1:sym_count
	index = symmetry(i)*2-1; % r-coordinate
	tmp = m(index,index);
	for j=1:msize
		m(j,index) = 0;
		m(index,j) = 0;
	end
	m(index,index) = tmp;
end

