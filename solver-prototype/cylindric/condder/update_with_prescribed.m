function nods = update_with_prescribed(nodes,boundary,dp)
	% Function for updating nodes with prescirbed displacements
	% and step dp
	nods = nodes;	
	boundary_count = size(boundary);
	boundary_count = boundary_count(1);

	if size(boundary,2) == 3
		boundary = convert_boundary(boundary);
	end	

	for i = 1:boundary_count
		if boundary(i,4) == 3 % all coordinates prescribed
			nods(boundary(i,1),1) = nods(boundary(i,1),1) + boundary(i,2);
			nods(boundary(i,1),2) = nods(boundary(i,1),2) + boundary(i,3);
		elseif boundary(i,4) == 1 % only r prescribed
			nods(boundary(i,1),1) = nods(boundary(i,1),1) + boundary(i,2);
		elseif boundary(i,4) == 2 % only z prescribed
			nods(boundary(i,1),2) = nods(boundary(i,1),2) + boundary(i,3);
		end
	end				
