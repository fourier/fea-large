function M = elasticity_matrix(F)
% function M = elasticity_matrix(F)
% Returns Elasticity matrix 4x4 for
% specified global g_model
% A5: lambda*I1(Epsilon)+2mu*Epsilon

	l = 0;
	m = 0;
	global g_model;
	global g_lambda;
	global g_mu;
	if strcmpi(g_model,'A5') % true if equal
		l = g_lambda;
		m = g_mu; 
	elseif strcmpi(g_model,'neohookean-compressible')
		J = det(F);
		l = g_lambda/J;
		m = (g_mu - g_lambda*log(J))/J;
	end			 
	M = [l+2*m, l, l, 0;
			 l,	l+2*m, l, 0;
			 l, l, l+2*m, 0;
			 0,	0,	0,	2*m;];

