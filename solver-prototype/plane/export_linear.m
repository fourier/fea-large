function export_linear(fname,elements,nodes,IDX,X)

fil = fopen(fname,'wt');
if fil == -1
	error('Unable to create file for output');
end

siz = size(elements);
elements_count = siz(1);
el_size = size(elements,2);

for i=1:elements_count
	if el_size == 3 % 3-noded triangle
		el = element(nodes,elements,i);
	else % 6-noded triangle
		el = element_full(nodes,elements,i);
	end

	v = zeros(el_size*2,1);
	v(1) = el(1,1); v(2) = el(1,2);
	v(3) = el(2,1); v(4) = el(2,2);
	v(5) = el(3,1); v(6) = el(3,2);
	
	A = elasticity_matrix();

	q = zeros(el_size*2,1);
	for j = 1:el_size*2
		q(j) = X(IDX(i,j));	
	end

	cntr = center_elem(el);

	d = v+q;
	Epsilon = epsilon_lin(el,q,cntr(1),cntr(2));
	Sigma = A*Epsilon;
	fprintf(fil,'%f,%f,%f,%f,%f,%f,',v(1),v(2),v(3),v(4),v(5),v(6));
	fprintf(fil,'%f,%f,%f,%f,%f,%f,0,',d(1),d(2),d(3),d(4),d(5),d(6));
	fprintf(fil,'%e,%e,0,%e\n',Sigma(1),Sigma(2),Sigma(3));
end

fclose(fil);

