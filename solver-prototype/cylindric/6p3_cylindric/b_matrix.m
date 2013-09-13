function B = b_matrix(el,r,z)
% Function for calculating B matrix by given element el 
% and point (r,z)

	B = zeros(4,12);
	L = l_matrix(el,r,z);

	% L(1,i) means L^i_{11}
	% L(2,i) means L^i_{22}
	% L(3,i) means L^i_{33}
	% L(4,i) means L^i_{12}
	% L(5,i) means L^i_{21}

	for i = 1:12
		B(1,i) = L(1,i);
		B(2,i) = L(2,i);
		B(3,i) = L(3,i);
		B(4,i) = 0.5*( L(4,i)+L(5,i) );
	end


