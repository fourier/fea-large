function el = element_full(nodes,elements,i)
% function el = element_full(nodes,elements,i)
% Function returns array of element's coordinates
% in form ((r1,z1),(r2,z2),(r3,z3),(r4,z4),(r5,z5),(r6,z6))

	el = [nodes(elements(i,1),1),nodes(elements(i,1),2);
				nodes(elements(i,2),1),nodes(elements(i,2),2);
				nodes(elements(i,3),1),nodes(elements(i,3),2);
				nodes(elements(i,4),1),nodes(elements(i,4),2);
				nodes(elements(i,5),1),nodes(elements(i,5),2);
				nodes(elements(i,6),1),nodes(elements(i,6),2)];

	
