function K = construct_kconstr(el0,el)
% function K = construct_kconstr(el0,el)
% Function for constructing Kconstr matrix 
% by given element el0 and current element el
K = zeros(12,12);
K = integrate7initial(el0,el,@kconstr);

function kc = kconstr(el0,el,r0,z0,r,z)
	q = create_q_from_elements(el0,el);
  % deformation gradient
  F = graddef(el0,r0,z0,q);
	J = det(F);
	if J <= 0.0
		string = sprintf('construct_kconstr: warning, gradF = %f', J);
		logger(string); 
%		stop_calculation_because_of_gradient
	end
%	point_initial  = [r;z;0];
%	point_current = F*point_initial;
%	r = point_current(1);
%	z = point_current(2);

	B = b_matrix(el,r,z);
	
	A = elasticity_matrix(F);
	kc = B'*A*B;
