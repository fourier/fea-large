function el = element(nodes,elements,i)
% Function returns array of element's coordinates
% in form ((x1,y1),(x2,y2),(x3,y3))

	el = [nodes(elements(i,1),1),nodes(elements(i,1),2);
				nodes(elements(i,2),1),nodes(elements(i,2),2);
				nodes(elements(i,3),1),nodes(elements(i,3),2)];
	
