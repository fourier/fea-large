function F = graddef2(full_element0,full_element,r,z)
	% function F = graddef2(full_element0,full_element,r,z)
	% Function for calculation Gradient of Deformations F
	% in current element full_element in point (r,z) where both 
	% full_element and (r,z) are given in current configuraion	
	% using initial element full_element0
	% F = dx/dX	
	% F_{ij} = \delta_{ij}+\tilde{u}_{ji} = 
	% = \delta+{ij}+\sum\limits_{k=1}^{12}L^k_{ji}*q_k
	x = zeros(3,6);
	for j = 1:6
		x(1,j) = full_element0(j,1);
		x(2,j) = full_element0(j,2);
		x(3,j)   = 0;
	end
	
	dform = zeros(3,6);
	for j = 1:6
		dform(1,j) = dform_dr(full_element,j,r,z);
		dform(2,j) = dform_dz(full_element,j,r,z);
		dform(3,j) = form(full_element,j,r,z)/r;
	end	
	
	F = zeros(3,3);
		for i = 1 : 3
			for j = 1 : 3
				sum = 0;
				for k = 1 : 6
					sum = sum + dform(j,k)*x(i,k);
				end			
				F(i,j) = sum;
			end
		end
	F(3,3) = 1;
	F = inv(F);
