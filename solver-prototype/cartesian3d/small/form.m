function f = form(el,i,x,y,z)
% f = form(el,i,x,y,z)
% Function to calculate shape function of element el
% by given point (x,y,z)
% i - index of shape function (1..10)
	f = 0;
%	l1 = local(el,1,x,y,z);
%	l2 = local(el,2,x,y,z);
%	l3 = local(el,3,x,y,z);
%	l4 = local(el,3,x,y,z);

	if ( i < 5 )
		f = (2*local(el,i,x,y,z)-1)*local(el,i,x,y,z);
	elseif i == 5 
		f = 4*local(el,1,x,y,z)*local(el,2,x,y,z);
	elseif i == 6
		f = 4*local(el,2,x,y,z)*local(el,3,x,y,z);
	elseif i == 7
		f = 4*local(el,3,x,y,z)*local(el,1,x,y,z);
	elseif i == 8
		f = 4*local(el,1,x,y,z)*local(el,4,x,y,z);
	elseif i == 9
		f = 4*local(el,2,x,y,z)*local(el,4,x,y,z);
	else % i == 10
		f = 4*local(el,3,x,y,z)*local(el,4,x,y,z);
	end

