function B = b_matrix6(el,x,y)
% Function for calculating B matrix by given element el 
% and point (r,z)

	% initialize
	B = zeros(3,12);
	
	for i = 1:6
		B(1,i*2-1) = dform6_dx(el,i,x,y);
		B(2,i*2) = dform6_dy(el,i,x,y);
		B(3,i*2-1) = dform6_dy(el,i,x,y);
		B(3,i*2) = dform6_dx(el,i,x,y);
	end
