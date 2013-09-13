function [E,S] = stress_in_nodes(nodes,elements,IDX,index,X)
% function s = stress(nodes,elements,IDX,index,X)
q = create_q_from_u(IDX,X,index);
el = element(nodes,elements,index);

nodes_count = size(elements,2);

E = zeros(6,nodes_count);
S = zeros(6,nodes_count);


for i = 1:nodes_count
	r = 0; s = 0; t = 0;
	if i == 1
		r = 0; s = 0; t = 0;	    
	elseif i == 2
	    r = 1; s = 0; t = 0;
	elseif i == 3
	    r = 0; s = 1; t = 0;
	elseif i == 4
	    r = 0; s = 0; t = 1;
	elseif i == 5
	    r = 0.5; s = 0; t = 0;
	elseif i == 6
		r = 0.5; s = 0.5; t = 0;
	elseif i == 7
		r = 0; s = 0.5; t = 0;
	elseif i == 8
		r = 0; s = 0; t = 0.5;
	elseif i == 9
		r = 0.5; s = 0; t = 0.5;
	else % i == 10
		r = 0; s = 0.5; t = 0.5;
	end

	e = epsilon_lin(el,q,r,s,t);
    E(:,i) = e;
	A = elasticity_matrix();
	S(:,i) = A*e;
end
