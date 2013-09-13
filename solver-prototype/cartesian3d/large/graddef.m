function F = graddef(el0,r,s,t,q)
	% function F = graddef(el0,r,s,t,q)
	% Function for calculation Gradient of Deformations F
	% in current element el in point (r,s,t) using 
	% already found nodal displacements q
	% 
	% F_{ij} = \sum\limits_{a=1}^{n} x_{a,i}\Der{N_a}{X_j}

	
	elesize = size(el0,1);
	F = zeros(3,3);
	dNs = dforms(el0,r,s,t);
	for i = 1 : 3
		for j = 1 : 3
			sum = 0;
			for k = 1 : elesize
				sum = sum + (el0(k,i)+q((k-1)*3+i))*dNs(j,k);
			end			
			F(i,j) = sum;
		end
	end

