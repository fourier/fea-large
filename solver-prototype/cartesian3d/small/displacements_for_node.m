function [ux,uy,uz] = displacements_for_node(X,node_number)
	ux = X(node_number*3-2);
	uy = X(node_number*3-1);
	uz = X(node_number*3-0);
end
