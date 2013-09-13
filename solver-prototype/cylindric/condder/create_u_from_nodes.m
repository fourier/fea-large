function U = create_u_from_nodes(msize,nodes0,nodes,elements,IDX)
	% function u = create_u_from_nodes(msize,nodes0,nodes,elements,IDX)
	% u = nodes0 - nodes	
	elements_count = size(elements,1);
	nodes_count = size(elements,2);
	U = zeros(msize,1);
	for i = 1:elements_count
		for j = 1:nodes_count			
			U(IDX(i,j*2-1)) = nodes(elements(i,j),1) - nodes0(elements(i,j),1);
			U(IDX(i,j*2)) = nodes(elements(i,j),2) - nodes0(elements(i,j),2);
		end
	end	

