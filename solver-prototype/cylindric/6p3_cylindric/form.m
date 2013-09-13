function f = form(el,i,r,z)
% Function to calculate shape function of element el
% by given point (r,z)
% i - index of shape function (1..6)
	f = 0;
	
	if ( i < 4 )
		f = (2*local(el,i,r,z)-1)*local(el,i,r,z);
	elseif i == 4 
		f	= 4*local(el,1,r,z)*local(el,2,r,z);
	elseif i == 5
		f = 4*local(el,2,r,z)*local(el,3,r,z);
	else % i == 6
		f = 4*local(el,3,r,z)*local(el,1,r,z);
	end
