% Program for calculation uniaxial stresses from Neo-Hookean model
% described in "Nonlinear Continuum Mechanics for Finite Element Analysis"
% by J.Bonet and R.D.Wood, second edition
% page 163

% Neo-Hookean material parameters - mu and lambda
global mu;
global lambda;
mu = 100;
lambda = 100;

% auxulary global K1 for solving non-linear equation
global K1;

% material model for Neo-Hookean material is 
% T = (mu/J)*(B-E)+(lambda/J)*lnJ*E
% where B is a Finger stress tensor, B = F*F^T
% J = det(F)

function s = k2_from_sigma22(k2)
	global K1;
	global mu;
	global lambda;
	J = K1*k2*k2;
	s = (mu*(k2**2-1)+lambda*log(J))/J; 	
end

function t = Txx(k1)
	global lambda;
  global mu;
	global K1;
	K1 = k1; 
	[k2,info] = fsolve(@k2_from_sigma22,1); 
	J = k1*k2*k2;
	t = (mu*(k1**2-1)+lambda*log(J))/J;
%	fprintf('k1 = %f, k2 = %f, J = %f, T = %f\n',k1,k2,J,t);
end

k = 1;
step = 0.008333;
for i=1:121
    stress = Txx(k);
    fprintf('%f %f\n',(k-1)*100,stress);
    k = k+step;
end
