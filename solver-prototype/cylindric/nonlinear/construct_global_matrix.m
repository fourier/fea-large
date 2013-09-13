function [M,P] = construct_global_matrix(msize,nodes,elements,IDX,X)
% Function for constructing global nonlinear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix,
% X - vector of previous solution
global g_compressible;

M = zeros(msize);
M = sparse(M);

P = zeros(msize,1);
siz = size(elements);
elements_count = siz(1);
y = zeros(3,1);

for i = 1:elements_count
	el = element(nodes,elements,i);
	S = s_elem(el);
	
	q = zeros(12,1);
	for j = 1:12
		q(j) = X(IDX(i,j));	
	end

	if not(g_compressible) % Incompressible case

		for j = 1:3
			y(j) = X(IDX(i,12+j));
		end
		MZE = construct_mze(el,q,y);
		% distribute in global matrix
		for k = 1:15
			for l = 1:15
				M(IDX(i,k),IDX(i,l))=M(IDX(i,k),IDX(i,l))+MZE(k,l);
			end
		end

		% distribution of boundary conditions in form of stresses
		h = 2*pi*S*integrate13(el,@construct_h,q);
		for k = 1:3
			P(IDX(i,12+k))=P(IDX(i,12+k))+h(k);		
		end
	else % Compressible case
		MZE = construct_mze(el,q,y);
		% distribute in global matrix
		for k = 1:12
			for l = 1:12
				M(IDX(i,k),IDX(i,l))=M(IDX(i,k),IDX(i,l))+MZE(k,l);
			end
		end
	end

%	fprintf('Stiffnes element %d from %d done...\n',i,elements_count);
end % finished loop through elements

