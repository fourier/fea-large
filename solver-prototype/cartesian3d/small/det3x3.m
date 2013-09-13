function d = det3x3(m)
% function d = det3x3(m)
% function for calculating determinant of matrix 3x3
  d = m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2)) - ...
  		m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1)) + ...
	  	m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1));
  
