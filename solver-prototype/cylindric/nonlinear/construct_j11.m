function j11 = construct_j11(el,q,y)
global g_compressible;

S = s_elem(el);
h = 0.0001;

K = 2*pi*S*integrate7(el,@construct_k,q);

if not(g_compressible)
	N = 2*pi*S*integrate7(el,@construct_n,q);
	f0 = K*q - N*y;
else
	f0 = K*q;
end


% Simple derivative

for j = 1:12
	qj = q;
	qj(j) = qj(j) + h;
	Kj = 2*pi*S*integrate7(el,@construct_k,qj);
	if not(g_compressible)
		Nj = 2*pi*S*integrate7(el,@construct_n,qj);
		fi = Kj*qj-Nj*y;
	else
		fi = Kj*qj;
	end
	j11(:,j) = (fi - f0)/h;
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
% Here implemented central derivative of order o(h^2)
%
for j = 1:12
	qj = q;
	qj(j) = q(j) + h;
	Kj = 2*pi*S*integrate7(el,@construct_k,qj);
	Nj = 2*pi*S*integrate7(el,@construct_n,qj);
	fi1 = Kj*qj-Nj*y;
	qj(j) = q(j) - h;
	Kj = 2*pi*S*integrate7(el,@construct_k,qj);
	Nj = 2*pi*S*integrate7(el,@construct_n,qj);
	fi2 = Kj*qj-Nj*y;
	j11(:,j) =  (fi1 - fi2)/(2*h);
end

