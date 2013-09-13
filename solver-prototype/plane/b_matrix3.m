function B = b_matrix3(el,x,y)
% Function for calculating B matrix by given element el 
% and point (r,z)

	% initialize
	B = zeros(3,6);
	% get size of element
	A = s_elem(el);
	
	% get coordintates of nodes
	x1 = el(1,1);
	x2 = el(2,1);
	x3 = el(3,1);
	y1 = el(1,2);
	y2 = el(2,2);
	y3 = el(3,2);

	% get "derivatives"
	y23 = y2-y3;
	y31 = y3-y1;
	y12 = y1-y2;
	x32 = x3-x2;
	x13 = x1-x3;
	x21 = x2-x1;

	% fill B matrix ...
	% first row of B matrix
	B(1,1) = y23; B(1,3) = y31; B(1,5) = y12;
	% second row of B matrix
	B(2,2) = x32; B(2,4) = x13; B(2,6) = x21;
	% third row of B matrix
	B(3,1) = x32; B(3,2) = y23; B(3,3) = x13; B(3,4) = y31; B(3,5) = x21; B(3,6) = y12;

	% ... and divide to 2*A
	B = 0.5*B/A;
