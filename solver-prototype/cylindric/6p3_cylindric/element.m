function el = element(nodes,elements,i)
% Function returns array of element's coordinates
% in form ((r1,z1),(r2,z2),(r3,z3))

	el = [nodes(elements(i,1),1),nodes(elements(i,1),2);
				nodes(elements(i,2),1),nodes(elements(i,2),2);
				nodes(elements(i,3),1),nodes(elements(i,3),2)];
	
