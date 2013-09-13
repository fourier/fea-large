function IDX = index_matrix(elements)
% function IDX = index_matrix(elements)
% calculates index matrix by given elements array

% Number of elements
elements_count = size(elements,1);
% Number of nodes per element
nodes_per_element = size(elements,2);
dof = 3;

IDX = zeros(elements_count,nodes_per_element*dof);

for k = 1:elements_count % number of elements
	for l = 1:nodes_per_element % number of nodes per element
    IDX(k,l*dof-2) = elements(k,l)*dof-2;
		IDX(k,l*dof-1) = elements(k,l)*dof-1;
		IDX(k,l*dof)   = elements(k,l)*dof;
	end
end
