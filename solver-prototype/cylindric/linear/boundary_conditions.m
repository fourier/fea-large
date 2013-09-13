function [m,p] = boundary_conditions(M,P,boundary,symmetry,dp)
% function [m,p] = boundary_conditions(M,P,boundary,symmetry,dp)
% Apply boundary conditions to the stiffness matrix M and 
% right-side vector P
% Parameters:
% M - global stiffness matrix
% P - global vector of external loads
% boundary - matrix of boundary conditions
% symmetry - matrix of nodes under symmentry condition(r- fixed)
% dp - additive coefficient for boundary conditions,
% used in iterative methods for solution of nonlinear tasks

m = M;
p = P;

boundary_count = size(boundary,1);
msize = length(m);

if size(boundary,2) == 3
	boundary = convert_boundary(boundary);
end	

% apply r-coordinate conditions
for i = 1:boundary_count
	% if it is "fake" condtion(0) or z-condition only - skip
	if boundary(i,4) == 0 || boundary(i,4) == 2
		continue
	end	
	index = boundary(i,1)*2-1; % r-coordinate
	tmp = m(index,index);
	condition = boundary(i,2);
	if ~(condition == 0)
		condition = condition + dp;
	end
	for j = 1:msize
		
		p(j) = p(j) - m(index,j)*condition;
		m(index,j) = 0;
		m(j,index) = 0;
	end
	m(index,index) = tmp;
	p(index,1) = tmp*condition;
end
% apply z-coordinate conditions
for i = 1:boundary_count
	% if it is "fake" condtion(0) or r-condition only - skip
	if boundary(i,4) == 0 || boundary(i,4) == 1
		continue
	end	
	index = boundary(i,1)*2; % z-coordinate
	tmp = m(index,index);
	condition = boundary(i,3);
	if ~(condition == 0)
		condition = condition + dp;
	end
	for j = 1:msize
		p(j) = p(j) - m(index,j)*condition;
		m(index,j) = 0;
		m(j,index) = 0;
	end
	m(index,index) = tmp;
	p(index,1) = tmp*condition;
end
% apply symmetry conditions
sym_count = size(symmetry,2);
for i=1:sym_count
	index = symmetry(i)*2-1; % r-coordinate
	tmp = m(index,index);
	for j=1:msize
		m(j,index) = 0;
		m(index,j) = 0;
	end
	m(index,index) = tmp;
	p(index,1) = 0;
end

