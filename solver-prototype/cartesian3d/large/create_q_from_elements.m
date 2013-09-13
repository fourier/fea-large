function q = create_q_from_elements(el0,el)
	% function q = create_q_from_elements(el0,el)
	% Creates displacements vector from elements el0 and el
	elesize = size(el,1);
	q = zeros(elesize,1);
	for i = 1:elesize
		q(i*3-2) = el(i,1) - el0(i,1);
		q(i*3-1) = el(i,2) - el0(i,2);
		q(i*3) = el(i,3) - el0(i,3);
	end	
