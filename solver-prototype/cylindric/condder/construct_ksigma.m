function ks = construct_ksigma(el0,el)
	% Function for constuction geometry stiffness matrix
	% of element el
	% el - element in current configuration
	% Ksigma = int ( \gradN * Sigma * \gradN )
	ks = zeros(12,12);
	ks = integrate7initial(el0,el,@ksigma);

function ks = ksigma(el0,el,r0,z0,r,z)
	ks = zeros(12,12);
	q = create_q_from_elements(el0,el);
  % deformation gradient
  F = graddef(el0,r0,z0,q);
	J = det(F);
	if J <= 0.0
		string = sprintf('construct_ksigma: warning, gradF = %f',J);
		logger(string); 
		%stop_calculation_because_of_gradient
	end	

%	point_initial  = [r;z;0];
%	point_current = F*point_initial;
%	r = point_current(1);
%	z = point_current(2);
	Sigma = sigma(F);
	E = eye(2,2);
	for k = 1:6
		for l = 1:6
			Bk = b_vector(el,k,r,z);
			Bl = b_vector(el,l,r,z);
			Hkl = Bk'*Sigma*Bl*E;
			ks(k*2-1:k*2,l*2-1:l*2) = Hkl;
		end
	end	
