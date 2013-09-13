function j21 = construct_j21(el,q,y)

S = s_elem(el);
h = 0.0001;

Bp = 2*pi*S*integrate7(el,@construct_bp,q);
H = 2*pi*S*integrate7(el,@construct_h,q);

f0 = -Bp*q - H;

% Simple derivative
for j = 1:12
	qj = q;
	qj(j) = qj(j) + h;
	Bpj = 2*pi*S*integrate7(el,@construct_bp,qj);
	Hj = 2*pi*S*integrate7(el,@construct_h,qj);
	fi = -Bpj*qj-Hj;
	j21(:,j) = (fi - f0)/h;
end
return
% Forward derivative of order o(h^2)
% 
% df    -f(x+2*h)+4*f(x+h)-3*f(x)
% -- =  -------------------------
% dx               2*h
%
% Backward derivative of order o(h^2)
%
% df   3*f(x)-4*f(x-h)+f(x-2*h)
% -- = ------------------------
% dx               2*h
% 
% Central derivative of order o(h^2):
%
% df    f(x+h)-f(x-h)
% -- =  -------------
% dx         2*h
%
% Here implented central derivative of order o(h^2)
%
for j = 1:12
	qj = q;
	qj(j) = q(j) + h;
	Bpj = 2*pi*S*integrate7(el,@construct_bp,qj);
	Hj = 2*pi*S*integrate7(el,@construct_h,qj);
	fi1 = -Bpj*qj-Hj;
	qj(j) = q(j) - h;
	Bpj = 2*pi*S*integrate7(el,@construct_bp,qj);
	Hj = 2*pi*S*integrate7(el,@construct_h,qj);
	fi2 = -Bpj*qj-Hj;
	j21(:,j) =  (fi1 - fi2)/(2*h);
end
