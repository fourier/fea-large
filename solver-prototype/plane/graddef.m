function F = graddef(el,x,y,q)
	% function F = graddef(el,r,z,q)
	% Function for calculation Gradient of Deformations F
	% in current element el in point (r,z) using 
	% already found nodal displacements q
	% 
	% F_{ij} = \delta_{ij}+\tilde{u}_{ji} = 
	% = \delta_{ij}+\sum\limits_{k=1}^{12}dNj/dxi*q_k
	F = eye(3,3);
	for i = 1 : 2
		for j = 1 : 2
			sum = 0;
			for k = 1 : 6
				if j == 1 % x
					sum = sum + dform6_dx(el,k,x,y)*q(2*(k-1)+i);
				else % j == 2 , y	
					sum = sum + dform6_dy(el,k,x,y)*q(2*(k-1)+i);
				end
			end			
			F(i,j) = F(i,j) + sum;
		end
	end

