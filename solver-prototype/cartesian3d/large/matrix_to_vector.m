function e = matrix_to_vector(E)
% e = matrix_to_vector(E)

	e = zeros(6,1);
	e(1) = E(1,1);
	e(2) = E(2,2);
	e(3) = E(3,3);
	e(4) = E(1,2);
	e(5) = E(2,3);
	e(6) = E(1,3);

