function L = l_matrix(el,r,z)
	% Function returns all values for L components
	% L = matrix(5,12), where column is
	% (L_11, L_22, L_33, L_12, L21)^T

	% L^i_{11} = dF_{1i}/dr
	% L^i_{22} = dF_{2i}/dz
	% L^i_{33} = F_{1i}/r
	% L^i_{12} = dF_{2i}/dr
	% L^i_{21} = dF_{1i}/dz
	% F_{ij} is a components of matrix of shape functions:
	%
	% F1 00	F2 00 F3 00 F4 00 F5 00 F6 00
	% 00 F1 00 F3 00 F3 00 F4 00 F5 00 F6
	
	L = zeros(5,12); % matrix of L components L11,L22,L33,L12,L21
	
	for i = 1:12	
		if rem(i,2) == 1
			j = ceil(i/2);
			L(1,i) = dform_dr(el,j,r,z);
			L(2,i) = 0;
			L(3,i) = form(el,j,r,z)/r;
			L(4,i) = 0;
			L(5,i) = dform_dz(el,j,r,z);
		else 
			j = floor(i/2);
			L(1,i) = 0;
			L(2,i) = dform_dz(el,j,r,z);
			L(3,i) = 0;
			L(4,i) = dform_dr(el,j,r,z);
			L(5,i) = 0;
		end
	end

