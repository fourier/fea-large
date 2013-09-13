function b = convert_boundary(boundary)
% function b = convert_boundary(boundary)
% converts boundary condtions to a new format
	b = boundary;
	boundary_count = size(boundary,1);
	if size(boundary,2) == 3
		% old boundary conditions, make them all of type 3
		b = zeros(boundary_count,4);
		b(1:boundary_count,1:3) = boundary;
		b(1:boundary_count,4) = 3;
	end	
	
