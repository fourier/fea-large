function el = element(nodes,elements,i)
% function el = element_full(nodes,elements,i)
% Function returns array of element's coordinates
% in form ((x1,y1,z1),...,(x10,y10,z10))

	el = [nodes(elements(i,1),1),nodes(elements(i,1),2),nodes(elements(i,1),3);
		  nodes(elements(i,2),1),nodes(elements(i,2),2),nodes(elements(i,2),3);
		  nodes(elements(i,3),1),nodes(elements(i,3),2),nodes(elements(i,3),3);
		  nodes(elements(i,4),1),nodes(elements(i,4),2),nodes(elements(i,4),3);
		  nodes(elements(i,5),1),nodes(elements(i,5),2),nodes(elements(i,5),3);
		  nodes(elements(i,6),1),nodes(elements(i,6),2),nodes(elements(i,6),3);
		  nodes(elements(i,7),1),nodes(elements(i,7),2),nodes(elements(i,7),3);
		  nodes(elements(i,8),1),nodes(elements(i,8),2),nodes(elements(i,8),3);
		  nodes(elements(i,9),1),nodes(elements(i,9),2),nodes(elements(i,9),3);
		  nodes(elements(i,10),1),nodes(elements(i,10),2),nodes(elements(i,10),3);];

	
