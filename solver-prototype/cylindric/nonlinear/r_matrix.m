function R = r_matrix(epsilon)
% Function for calculating R-matrix by given deformations vector

	R = zeros(4,4);
	R(1,2) = epsilon(3);
	R(2,1) = epsilon(3);
	R(3,1) = epsilon(2);
	R(3,4) = -epsilon(4);
	R(4,4) = -epsilon(3);
