function [k2_an,s1_an,k2_bn,s2_bn] = uniaxial(k1,n,lambda,mu,betta)
% function [k2_an,s1_an,k2_bn,s2_bn] = uniaxial(k1,n,lambda,mu,betta)
% Program for calculation stresses and strains for uniaxial
% case, models A1-A5 && B1-B5
	k2_an = k2_An(k1,lambda,mu,n)^(1./(n-3.));
	k2_bn = k2_Bn(k1,betta,mu,n)^(1./(n-3.));
	s1_an = sigma_an(k1,lambda,mu,n);
	s2_bn = sigma_bn(k1,betta,mu,n);	


function s = sigma_an(k1,l,m,n)
	% get k1^(n-3) which is used in equations
	kk1 = k1^(n-3);
	% obtain k2^(n-3) to substitute to the final equation
	kk2 = k2_An(k1,l,m,n);
	% obtain k2 from k2^(n-3)
	k2 = kk2^(1./(n-3.));
	% get sigma11
	% Sigma_11 = F_11^(n-3)*T_11/detF
	s = (k1^(n-3-1)/k2^2)*((l+2*m)*kk1+2*l*kk2-(3*l+2*m))/(n-3);
	% for small deformations F_11=1,detF=1
	% so the next formula is used
	%	s = ((l+2*m)*kk1+2*l*kk2-(3*l+2*m))/(n-3);

function s = sigma_bn(k1,b,m,n)
	kk1 = k1^(n-3);
	k2_pow_n_minus_3 = k2_Bn(k1,b,m,n);
	kk2 = k2_pow_n_minus_3;
	s = m*(n-3)*kk1*(1+b+(1-b)*(kk1+2*kk2)/(n-3)-(1-b)*kk1/(n-3));
	
function k2 = k2_An(k1,l,m,n)
	% function k2 = k2_An(k1,l,m,n)
	% returns k2^(n-3) for model An
	kk1 = k1^(n-3);
	k2 = (3*l+2*m-l*kk1)/(2.*l+2.*m);

function k2 = k2_Bn(k1,b,m,n)
	% function k2 = k2_Bn(k1,b,m,n)
	% k1 == k1^(n-3)
	% returns k2^(n-3) for model Bn
	kk1 = k1^(n-3);
	k2 = ((1+b)*(n-3)+(1-b)*kk1)/(b-1);
	
	
	
