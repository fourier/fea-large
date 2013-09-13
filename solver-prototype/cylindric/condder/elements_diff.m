function q = elements_diff(elements,nodes0,nodes,i)
	% Function elements_diff(elements,nodes0,nodes,i)
	% returns vector q in format of displacements with
	% difference between element i from nodes0 and nodes
	% nodes arrays
	q = zeros(12,1);
	for j = 1:6
		q(j*2-1) = nodes(elements(i,j),1)-nodes0(elements(i,j),1);
		q(j*2)   = nodes(elements(i,j),2)-nodes0(elements(i,j),2);
	end				

