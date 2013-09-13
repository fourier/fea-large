function q = create_q_from_elements(el0,el)
	% function q = create_q_from_elements(el0,el)
	% Creates displacements vector from elements el0 and el
	q = zeros(12,1);
	for i = 1:size(el0,1);
		q(i*2-1) = el(i,1) - el0(i,1);
		q(i*2) = el(i,2) - el0(i,2);
	end	
