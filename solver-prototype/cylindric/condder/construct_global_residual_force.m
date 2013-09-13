function R = construct_global_residual_force(msize,nodes0,nodes,elements,IDX)  
	% R = construct_global_residual_force(msize,nodes,elements,IDX)  
	% function for constructing global nodal residual force vector R
	% R = T - F
	% where T - nodal internal forces, F - nodal external forces
	% when no external forces, R = T

	R = zeros(msize,1);
	elements_count = size(elements,1);

	for i = 1:elements_count
		el = element_full(nodes,elements,i);
		el0 = element_full(nodes0,elements,i);
		S = s_elem(el);
		T = 2*pi*S*construct_t(el0,el);
		% distribute in global vector
		for k = 1:12
			R(IDX(i,k))=R(IDX(i,k))+T(k);			
		end	
	end
