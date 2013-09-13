function Ke = construct_mze(el,q,y)
% Constructs local stiffness matrix by given element el and solutions q,y
	global g_compressible;

	S = s_elem(el);

	K = 2*pi*S*integrate13(el,@construct_k,q);
	if not(g_compressible)
		N = 2*pi*S*integrate13(el,@construct_n,q);
		Bp = 2*pi*S*integrate13(el,@construct_bp,q);

		Ke = zeros(15,15);
		Ke(1:12,1:12) = K;
		Ke(1:12,13:15) = -N;
		Ke(13:15,1:12) = -Bp;
	else
	    Ke = zeros(12,12);
		Ke = K;
	end

