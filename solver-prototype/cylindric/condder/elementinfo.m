function [F,E,Elin,Sigma,q] = elementinfo(nodes,elements,IDX,X,index)
	% function [F,E,Elin,Sigma,q] = elementinfo(nodes,elements,IDX,X,index)
	% Returns information about element with index 'index' in center
	% of element: 
	% F - Deformation gradient in center
	% E - Right Cauchy-Green strain tensor
	% Elin - tesnor of small strains
	% Sigma - Cauchy stress tensor
	% q - displacements in all nodes (1-12)
	el = element(nodes,elements, index);
	q = create_q_from_u(IDX,X,index);
	cntr = center_elem(el);
	F = graddef(el,cntr(1),cntr(2),q);
	Sigma = sigma(F);	
	E = 0.5*(F'*F-eye(3));
	Elin = epsilon_lin(el,q,cntr(1),cntr(2));
