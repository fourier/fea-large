function [M,P] = construct_global_matrix(msize,nodes,elements,IDX)
% Function for constructing global linear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix,

% initialize global matrix of stiffness
M = zeros(msize);
M = sparse(M);

% global vector of external forces
P = zeros(msize,1);


elements_count = size(elements,1);
siz = size(elements,2)*2;

global indd;
indd = 1;
% loop through elements
for i = 1:elements_count

	if size(elements,2) == 3 % 3-noded triangle
		el = element(nodes,elements,i);
	else % 6-noded triangle
		el = element_full(nodes,elements,i);
	end
	S = s_elem(el);
	cntr = center_elem(el);

	% integrate in central point
%	K = S*construct_k(el,cntr(1),cntr(2));
%	K = S*integrate7(el,@construct_k);
	K = S*integrate13(el,@construct_k);
	fname = sprintf('Dump/K%d.txt',i);
	dumpmatrix(fname,K);
	% distribute in global matrix
	for k = 1:siz
		for l = 1:siz
			M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + K(k,l);
		end
	end
end


