function S = sigma(F)
	% function S = sigma(F)
	% Function for calculation Cauchy stress
	% F - deformation gradient matrix
	% return S 3x3 matrix
	S = zeros(3,3);
	global g_model;
	if strcmpi(g_model,'A5')
		S = sigma_a5(F);
	elseif strcmpi(g_model,'neohookean-compressible')
		S = sigma_neohookean_compressible(F);	
	end	
		
