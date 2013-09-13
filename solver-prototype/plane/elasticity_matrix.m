function M = elasticity_matrix()
% function M = elasticity_matrix(F)
% Returns Elasticity matrix 4x4 for
% specified global g_model
%	E = 1e9;
%	nu = 0.3;

l = 100;
m = 100;
E = m*(3*l+2*m)/(l+m);
nu = 0.5*l/(l+m);
%	M(1,1) = E/(1-nu^2);
%	M(2,2) = M(1,1);
%	M(3,3) = 0.5*E/(1+nu);
%	M(1,2) = nu*M(1,1);
%	M(2,1) = M(1,2);

M = [l+2*m, l, 0;
	   l, l+2*m, 0;
		 0, 0, m];
