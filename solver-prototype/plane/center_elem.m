function c = center_elem(el)
% Function for calculating center of triangle element

	c = [1/3.0*(el(1,1)+el(2,1)+el(3,1)),1/3.0*(el(1,2)+el(2,2)+el(3,2))];

