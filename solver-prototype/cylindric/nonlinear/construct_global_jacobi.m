function J = construct_global_jacobi(msize,nodes,elements,IDX,X)
% Function for constructing global nonlinear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix,
% X - vector of previous solution
global g_compressible;

J = zeros(msize);
J = sparse(J);

siz = size(elements);
elements_count = siz(1);

for i = 1:elements_count
	el = element(nodes,elements,i);
	S = s_elem(el);
	
	q = zeros(12,1);
	for j = 1:12
		q(j) = X(IDX(i,j));	
	end
	y = zeros(3,1);
	
	if not(g_compressible)

		for j = 1:3
			y(j) = X(IDX(i,12+j));
		end

		N = 2*pi*S*integrate13(el,@construct_n,q);

		Ke = zeros(15,15);
		Ke(1:12,1:12) = construct_j11(el,q,y);
		Ke(1:12,13:15) = -N;
		Ke(13:15,1:12) = construct_j21(el,q,y);

		% distribute in global matrix
		for k = 1:15
			for l = 1:15
				J(IDX(i,k),IDX(i,l))=J(IDX(i,k),IDX(i,l))+Ke(k,l);
			end
		end
	else % Compressible case
		Ke = zeros(12,12);
		Ke = construct_j11(el,q,y);
		for k = 1:12
			for l = 1:12
				J(IDX(i,k),IDX(i,l))=J(IDX(i,k),IDX(i,l))+Ke(k,l);
			end
		end
	end

%	fprintf('Jacobi element %d from %d done...\n',i,elements_count);
end % finished loop through elements

