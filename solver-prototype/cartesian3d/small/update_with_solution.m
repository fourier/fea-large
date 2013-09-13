function nods = update_with_solution(nodes,X)
	% Function for updating nodes with prescirbed displacements
	nodes_count = size(nodes,1);
		
	nods = zeros(nodes_count,3);	
	
	for i = 1:nodes_count
		nods(i,1) = nodes(i,1) + X(i*3-2);
		nods(i,2) = nodes(i,2) + X(i*3-1);
		nods(i,3) = nodes(i,3) + X(i*3);
	end				

