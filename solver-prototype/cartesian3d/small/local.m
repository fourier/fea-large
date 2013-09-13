function l = local(el,i,x,y,z)
% l = local(el,i,x,y,z)
% i = 1..4
% Function for calculating local tethraedra coordinates L1,L2,L3,L4
% by given index i and element el for point (x,y,z)

% V1 = vol(P234) V2 = vol(P143)
% V3 = vol(P124) V4 = vol(P213)

	l = 0;
	p = [x,y,z];	
	p1 = el(1,:);
	p2 = el(2,:);
	p3 = el(3,:);
	p4 = el(4,:);

	V = volume(el);
		
	if ( i == 1 )
		V1 = volume3(p,p2,p3,p4);
		l = V1/V;
	elseif i == 2 
		V2 = volume3(p,p1,p4,p3);
		l = V2/V;
	elseif i == 3
		V3 = volume3(p,p1,p2,p4);
		l = V3/V;
	else % i == 4
		V4 = volume3(p,p2,p1,p3);
		l = V4/V;
	end
	