function F = graddef1(full_element0,full_element,r,z)
	% function F = graddef1(full_element0,full_element,r,z)
	% Function for calculation Gradient of Deformations F
	% in current element el in point (r,z) where both el
	% and (r,z) are given in current configuraion	
	% F = dx/dX	
	% F_{ij} = \delta_{ij}+\tilde{u}_{ji} = 
	% = \delta_{ij}+\sum\limits_{k=1}^{12}L^k_{ji}*q_k
	x0 = zeros(12,1);
	for j = 1:6
		x0(j*2-1) = full_element0(j,1);
		x0(j*2) = full_element0(j,2);
	end
	L = l_matrix(full_element,r,z);
	Finv = zeros(3,3);
		for i = 1 : 3
			for j = 1 : 3
				sum = 0;
				for k = 1 : 12
					if i == j  			
						sum = sum + L(i,k)*x0(k);
					elseif i == 1 && j == 2
						sum = sum + L(5,k)*x0(k);
					elseif i == 2 && j == 1
						sum = sum + L(4,k)*x0(k);
					end
				end			
				Finv(i,j) = sum;
			end
		end
	F = inv(Finv);
