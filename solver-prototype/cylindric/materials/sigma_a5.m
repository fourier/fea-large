function S = sigma_a5(F)
	%	function S = sigma_a5(F)
	% Function for calculation Cauchy stress
	% F - deformation gradient matrix
	% return S 3x3 matrix
	global g_lambda;
	global g_mu;
	G = F'*F;
	E = eye(3,3);
	C = 0.5*(G-E);
	% Sigma = lambda*I1(C)+2mu*C
	S = g_lambda*invariant_1st(C)*E+2*g_mu*C;
