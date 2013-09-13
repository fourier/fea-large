function T = construct_t(el0,el)
	% function T = construct_t(el0,el)
	% Function for constuction of internal force vector T
	% el - element in current configuration
	% q - displacements from initial configuration(for deformation gradient)
	% T = int(B^T * Sigma)
	T = zeros(12,1);
	T = integrate7initial(el0,el,@internal_forces);
function r = internal_forces(el0,el,r0,z0,r,z)
	% derivatives of shape functions

	q = create_q_from_elements(el0,el);
	% deformation gradient 
	F = graddef(el0,r0,z0,q);
	J = det(F);
	if J <= 0.0
		string = sprintf('construct_t: warning, gradF = %f', J);
		logger(string); 
%		stop_calculation_because_of_gradient
	end

%	point_initial  = [r;z;0];
%	point_current = F*point_initial;
%	r = point_current(1);
%	z = point_current(2);

	B = b_matrix(el,r,z);

	% Cauchy stress tensor
	S = sigma(F);
	% Vector of stresses
	s = zeros(4,1);
	s(1) = S(1,1);
	s(2) = S(2,2);
	s(3) = S(3,3);
	s(4) = S(1,2);
	r = B'*s;
