function l = local(el,i,r,z)
% Function for calculating local triangle coordinates L1,L2,L3
% by given index i and element l for point (r,z)
	l = 0;
	
	a = zeros(1,3);
	a(1) = el(2,1)*el(3,2)-el(3,1)*el(2,2);
	a(2) = el(3,1)*el(1,2)-el(1,1)*el(3,2);
	a(3) = el(1,1)*el(2,2)-el(2,1)*el(1,2);

	b = zeros(1,3);
	b(1) = el(2,2)-el(3,2);
	b(2) = el(3,2)-el(1,2);
	b(3) = el(1,2)-el(2,2);

	c = zeros(1,3);
	c(1) = el(3,1)-el(2,1);
	c(2) = el(1,1)-el(3,1);
	c(3) = el(2,1)-el(1,1);

	S = s_elem(el);
	l = (1.0/(2.0*S))*(a(i)+b(i)*r+c(i)*z);

