function r = invariant_1st(A)
	siz = size(A);
	r = 0;
	for i = 1:siz(1)
		r = r + A(i,i);
	end	
