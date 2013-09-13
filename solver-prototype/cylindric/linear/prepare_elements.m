% prepare matrixes and data for solution
function [elements_count,corners_count,msize,IDX] = prepare_elements(nodes,elements)
global g_compressible;

%A = elasticity_matrix();
% let's count number of nodes, elements and corners
siz1 = size(nodes);
siz2 = size(elements);

% Number of nodes
nodes_count = siz1(1);
% Number of elements
elements_count = siz2(1);

corners = zeros(1,nodes_count*3);
l = 0;
for j = 1:elements_count
	for k = 1:3
		l = l+1;
		corners(l) = elements(j,k);
	end
end

% final corners count
zero_elements = 0;
for j = 1:nodes_count*3;
	if corners(j) == 0
		zero_elements = zero_elements + 1;
	end
end

corners_count = length(unique(corners)) ; % -1 because of 0
if zero_elements > 0
	corners_count = corners_count - 1;
end

% size of global matrix is number of degreese of freedom by
% displacements in all and hydrostatic pressure in corners of elements
if g_compressible
    % just all nodes in compressible case
	msize = nodes_count*2;
else
    % all nodes + corner nodes, because p has 1st order of approximation
	msize = nodes_count*2+corners_count; 
end

logger('Creating INDEX matrix...');
%	create an index matrix with degrees of freedom I(m,n), where
% m - number of element, n - number of degree of freedom(global)
% siz2(2)*2 - number of displacements(u,v) per element
% 3 - number of hydrostatic pressure values

index_rows = 0;
if g_compressible
	index_rows = siz2(2)*2; 
else 
	siz2(2)*2+3;
end

IDX = zeros(elements_count,index_rows);

% fill index matrix with degrees of freedom by displacements
for k = 1:elements_count % number of elements
	for l = 1:siz2(2) % number of nodes per element
		IDX(k,l*2-1) = elements(k,l)*2-1;
		IDX(k,l*2)   = elements(k,l)*2;
	end
end

if not(g_compressible)
	% fill index matrix with hydrostatic pressure d.o.f
	% starting from number of displacements d.o.f + 1
	corners_IDX = zeros(elements_count,3);
	p = 0;
	for i = 1:nodes_count % number of nodes
		found = 0;
		for k = 1:elements_count % number of elements
			if found > 0
				break;
			end
			for l = 1:3 % corner elements
				if elements(k,l) == i % found corner index
					if corners_IDX(k,l) == 0 % position is 0
						p = p + 1;
						corners_IDX(k,l) = p;
						for j = 1: elements_count
							for m = 1:3
								if elements(k,l) == elements(j,m)
									corners_IDX(j,m) = p;
								end
							end
						end
						found = 1;	
						break;
					end
				end
			end
		end
	end

	p_start = nodes_count*2;
	for i = 1:elements_count % number of elements
		for j = 1:3
			IDX(i,12+j) = corners_IDX(i,j) + p_start;
		end
	end
end

% fill index matrix with degrees of freedom by hydrostatic pressure

string = sprintf('Number of elements: %d',elements_count);
logger(string);
string = sprintf('Number of nodes: %d',nodes_count);
logger(string);
string = sprintf('Size of global stiffness matrix: %dx%d',msize,msize);
logger(string);