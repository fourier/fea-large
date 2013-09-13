function [M,P] = construct_global_matrix_lin(msize,nodes,elements,IDX)
% Function for constructing global linear stiffness matrix
% and global vector of loads
% Params:
% msize - size of global stiffness matrix and loads vector
% nodes - nodes matrix
% elements - elements matrix
% IDX - index matrix,

global g_compressible;

% initialize global matrix of stiffness
M = zeros(msize);
M = sparse(M);

% global vector of external forces
P = zeros(msize,1);
siz = size(elements);

elements_count = siz(1);

if not(g_compressible)
	% loop through elements
	for i = 1:elements_count
		el = element(nodes,elements,i);
		S = s_elem(el);
		

		% integrate numerically 
		K = 2*pi*S*integrate7(el,@construct_k_lin,[]);
		N = 2*pi*S*integrate7(el,@construct_n_lin,[]);
	
		% form local stiffness matrix
		Ke = zeros(15,15);
		Ke(1:12,1:12) = K;
		Ke(1:12,13:15) = -N;
		Ke(13:15,1:12) = -N';
	%	Alpha = 100;
	%	Ke(13:15,13:15) = -1/Alpha*eye(3);
	
		% distribute in global matrix
		for k = 1:15
			for l = 1:15
				M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + Ke(k,l);
			end
		end

		% distribution of boundary conditions in form of stresses
		% not done currently

	end % finished loop through elements
else % compressible case
	% loop through elements
	for i = 1:elements_count
		el = element(nodes,elements,i);
		S = s_elem(el);

		% integrate numerically 
		K = 2*pi*S*integrate7(el,@construct_k_lin,[]);
		for k = 1:12
			for l = 1:12
				M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + K(k,l);
			end
		end

		% distribution of boundary conditions in form of stresses
		% not done currently

	end % finished loop through elements
end


