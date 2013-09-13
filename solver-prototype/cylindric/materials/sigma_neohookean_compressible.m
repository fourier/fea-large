function S = sigma_neohookean_compressible(F)
	% function S = sigma_neohookean_compressible(F)
	% Function for calculation Cauchy stress
	% F - deformation gradient matrix
	% return S 3x3 matrix
	global g_lambda;
	global g_mu;
	J = det(F);
	b = F*F';
	E = eye(3,3);
	C = b-E;
	% Sigma = lambda/J*(ln(J))E+mu/J*(b-E)
	S = g_mu/J*C+g_lambda/J*log(J)*E;
