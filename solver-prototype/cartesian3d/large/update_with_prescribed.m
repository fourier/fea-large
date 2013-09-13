function nods = update_with_prescribed(nodes0,nodes,boundary,IT)
	% nods = update_with_prescribed(nodes0,nodes,boundary,IT)
	% Function for updating nodes with prescirbed displacements
	nods = nodes;	
	boundary_count = size(boundary);
	boundary_count = boundary_count(1);

	for i = 1:boundary_count
		for j = 1:3
			nods(boundary(i,1),j) = nodes0(boundary(i,1),j) + IT*boundary(i,1+j);
		end
	end				
