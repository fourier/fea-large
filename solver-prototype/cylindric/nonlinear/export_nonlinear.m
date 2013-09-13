function export_nonlinear(fname,elements,nodes,IDX,X)
global g_compressible;

A = elasticity_matrix();
I0 = ones(4,1);
I0(4) = 0;

fil = fopen(fname,'wt');

siz = size(elements);
elements_count = siz(1);

for i=1:elements_count
	el = element(nodes,elements,i);
	v = zeros(12,1);
	v(1) = el(1,1); v(2) = el(1,2);
	v(3) = el(2,1); v(4) = el(2,2);
	v(5) = el(3,1); v(6) = el(3,2);

	q = zeros(12,1);
	for j = 1:12
		q(j) = X(IDX(i,j));	
	end

		
	y = zeros(3,1);
	if not(g_compressible)
		for j = 1:3
			y(j) = X(IDX(i,12+j));
		end
	end
	
	cntr = center_elem(el);

	Fp = fp_matrix(el,cntr(1),cntr(2));

	d = v+q;
	Epsilon = epsilon_nonlin(el,q,cntr(1),cntr(2));
	Sigma = A*Epsilon - I0*Fp'*y; % will be just A*Epsilon in compressible case
%	Sigma = stabilize_sigma(el,Sigma,cntr(1),cntr(2)); 
	fprintf(fil,'%f,%f,%f,%f,%f,%f,',v(1),v(2),v(3),v(4),v(5),v(6));
	fprintf(fil,'%f,%f,%f,%f,%f,%f,0,',d(1),d(2),d(3),d(4),d(5),d(6));
	fprintf(fil,'%f,%f,%f,%f\n',Sigma(1),Sigma(2),Sigma(3),Sigma(4));
end

fclose(fil);

