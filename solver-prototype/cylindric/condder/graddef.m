function F = graddef(el,r,z,q)
	% function F = graddef(el,r,z,q)
	% Function for calculation Gradient of Deformations F
	% in current element el in point (r,z) using 
	% already found nodal displacements q
	% 
	% F_{ij} = \delta_{ij}+\tilde{u}_{ji} = 
	% = \delta_{ij}+\sum\limits_{k=1}^{12}L^k_{ji}*q_k
	L = l_matrix(el,r,z);
	F = zeros(3,3);
	for i = 1 : 3
		for j = 1 : 3
			sum = 0;
			for k = 1 : 12
				if i == j  			
					sum = sum + L(i,k)*q(k);
				elseif i == 1 && j == 2
					sum = sum + L(5,k)*q(k);
				elseif i == 2 && j == 1
					sum = sum + L(4,k)*q(k);
				end
			end			
			F(i,j) = kroneker_delta(i,j) + sum;
		end
	end

