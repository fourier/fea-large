function S = sigma(F)
%	function S = sigma(F)
    S = sigma_A5(F);
%    S = sigma_neohookean(F);
end
  
function S = sigma_A5(F)
	%	function S = sigma(F)
	mu = 100;
	lambda = 100;
	G = F'*F;
	E = eye(3,3);
	C = 0.5*(G-E);
	% Sigma = lambda*I1(C)+2mu*C
	S = (lambda*invariant_1st(C)*E+2*mu*C)/det(F);
  S = F*S*F';
end

function S = sigma_neohookean(F)
	mu = 100;
	lambda = 100;
	B = F*F';
	E = eye(3,3);
  J = det(F);
  S = mu/J*(B-E)+lambda/J*log(J)*E;
end