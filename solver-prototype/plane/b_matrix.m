function B = b_matrix(el,x,y)
% Function for calculating B matrix by given element el 
% and point (r,z)

if size(el,1) == 3 % 3-noded triangle
	B = b_matrix3(el,x,y);
else % 6-noded triangle
	B = b_matrix6(el,x,y);
end