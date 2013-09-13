% prepare matrixes and data for solution
function [elements_count,msize,IDX] = prepare_elements(nodes,elements)

%A = elasticity_matrix();
% let's count number of nodes, elements and corners
siz1 = size(nodes);
siz2 = size(elements);

% Number of nodes
nodes_count = siz1(1);
% Number of elements
elements_count = siz2(1);

msize = nodes_count*2;

index_rows = 0;
index_rows = siz2(2)*2; 

IDX = zeros(elements_count,index_rows);

% fill index matrix with degrees of freedom by displacements
for k = 1:elements_count % number of elements
	for l = 1:siz2(2) % number of nodes per element
		IDX(k,l*2-1) = elements(k,l)*2-1;
		IDX(k,l*2)   = elements(k,l)*2;
	end
end