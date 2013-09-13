function M = construct_global_matrix_cond(msize,nodes0,nodes,elements,IDX)
% M = construct_global_matrix_cond(msize,nodes0,nodes,elements,IDX)
% Function for constructing global nonlinear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix,

M = zeros(msize);
M = sparse(M);

siz = size(elements);
elements_count = siz(1);

for i = 1:elements_count
	el = element_full(nodes,elements,i);
	el0 = element_full(nodes0,elements,i);
	S = s_elem(el);
	MZE = zeros(12,12);
	Kc = 2*pi*S*construct_kconstr(el0,el);
	Ks = 2*pi*S*construct_ksigma(el0,el);
	MZE = Kc + Ks;
	% distribute in global matrix
	for k = 1:12
		for l = 1:12
			M(IDX(i,k),IDX(i,l))=M(IDX(i,k),IDX(i,l))+MZE(k,l);
		end
	end

end % finished loop through elements

