function nods = update_with_solution(nodes,elements,IDX,X)
	% Function for updating nodes with prescirbed displacements
	% and step dp
	nods = nodes;	
	elements_count = size(elements);
	msize = size(X,1);
	for i = 1:msize/2
		nods(i,1) = nods(i,1) + X(i*2-1);
		nods(i,2) = nods(i,2) + X(i*2);
	end				

	return 				
	for i = 1:elements_count(1)
		for j = 1:elements_count(2)
			nods(elements(i,j),1) = nods(elements(i,j),1) + X(IDX(i,j*2-1));
			nods(elements(i,j),2) = nods(elements(i,j),2) + X(IDX(i,j*2));		
		end			
	end	
