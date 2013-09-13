function b = b_vector(el,i,r,z)
	% function b = b_vector(el,i,r,z)
	% returns vector of derivatives of shape functions
	% for element el of node i for point (r,z)
	% B = [df_i/dr,df_i/dz,df_i/d\phi]
	
	b = [dform_dr(el,i,r,z); 
			 dform_dz(el,i,r,z);
			 form(el,i,r,z)/r];
