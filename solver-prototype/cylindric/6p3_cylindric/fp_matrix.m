function fp = fp_matrix(el,r,z)
% Function to calculate Fp vector of shape
% functions of hydrostatic pressure by given element el and point (r,z)

	fp = zeros(3,1);
	for j = 1:3
		fp(j)  = form_p(el,j,r,z);
	end
