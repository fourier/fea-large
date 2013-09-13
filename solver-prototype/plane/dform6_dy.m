function f = dform6_dy(el,i,x,y)
% Function for calculating derivative of i-th 
% shape function by y:
% dF_i/dy(x,y)
	
	f = 0;
	r = x;
	z = y;

	r1 = el(1,1);
	r2 = el(2,1);
	r3 = el(3,1);

	S = s_elem(el);

	switch i
		case 1
			f = 0.5*(r3-r2)*(4*local(el,1,r,z)-1)/S;
		case 2
			f = 0.5*(r1-r3)*(4*local(el,2,r,z)-1)/S;
		case 3
			f = 0.5*(r2-r1)*(4*local(el,3,r,z)-1)/S;
		case 4
			f = 2*((r3-r2)*local(el,2,r,z)+(r1-r3)*local(el,1,r,z))/S;
		case 5
			f = 2*((r1-r3)*local(el,3,r,z)+(r2-r1)*local(el,2,r,z))/S;
		case 6
			f = 2*((r2-r1)*local(el,1,r,z)+(r3-r2)*local(el,3,r,z))/S;
	end
	
