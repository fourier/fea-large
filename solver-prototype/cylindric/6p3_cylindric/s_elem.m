function s = s_elem(el)
% Function returns square size of element
	s = 0;
	A = ones(3);
	A(:,2:3) = el(1:3,:);
	s = det(A)/2.0;	
