function [u,el0,el,el1,q0,q] = temp_element_values(msize,...
	nodes0,nodes,elements,IDX,X0,index)
	% [u,el0,el,el1,q0,q] = temp_element_values(msize,nodes0,nodes,elements,IDX,X0,index)	
	% u - solution vector(from initial configuration)
	% el0 - initial element
	% el - element in actual configuration from elements for nodes 'nodes'
	% el1 - element in actual configuraion from linear solution X0
	% q0 - solution vector for element from linear solution
	% q - solution vector for element from condder solution 'nodes'

	u = create_u_from_nodes(msize,nodes0,nodes,elements,IDX);
	el0 = element_full(nodes0,elements,index);
	el = element_full(nodes,elements,index);
	q = create_q_from_u(IDX,u,index);
	q0 = create_q_from_u(IDX,X0,index);
	el1 = el0;
	for i = 1:6
		for j = 1:2
		el1(i,j) = el1(i,j)+q0((i-1)*2+j);
		end
	end


