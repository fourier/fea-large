function m = inv3x3(a)
% function m = inv3x3(a)
% calculate inverse matrix 3x3
  m = zeros(3,3);
	d = det3x3(a);
	% calculate components
	% first row
	m(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/d;
	m(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/d;
	m(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))/d;
	% second row
	m(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/d;
	m(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))/d;
	m(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/d;
	% third row
	m(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/d;
	m(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/d;
	m(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/d;
  

