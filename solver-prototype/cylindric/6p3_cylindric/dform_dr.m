function f = dform_dr(el,i,r,z)
% function f = dform_dr(el,i,r,z)
% Function for calculating derivative of i-th 
% shape function by r:
% dF_i/dr(r,z)
	
	f = 0;

	z1 = el(1,2);
	z2 = el(2,2);
	z3 = el(3,2);

	S = s_elem(el);
	
	switch i
		case 1
			f = 0.5*(z2-z3)*(4*local(el,1,r,z)-1)/S;
		case 2
			f = 0.5*(z3-z1)*(4*local(el,2,r,z)-1)/S;
		case 3
			f = 0.5*(z1-z2)*(4*local(el,3,r,z)-1)/S;
		case 4
			f = 2*((z2-z3)*local(el,2,r,z)+(z3-z1)*local(el,1,r,z))/S;
		case 5
			f = 2*((z3-z1)*local(el,3,r,z)+(z1-z2)*local(el,2,r,z))/S;
		case 6
			f = 2*((z1-z2)*local(el,1,r,z)+(z2-z3)*local(el,3,r,z))/S;
	end
	
