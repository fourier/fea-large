function B = b1_matrix(el,q,r,z)
% Function for calculating \hat{B} matrix
% by given element and point(r,z)

	B = zeros(4,12);

	% i-column of \hat{B} matrix is
  % row1 = \sum\limits^{12}_{j=1}H^{ji}_{11}q_j
	% row2 = \sum\limits^{12}_{j=1}H^{ji}_{22}q_j
	% row3 = \sum\limits^{12}_{j=1}H^{ji}_{33}q_j
	% row4 = \sum\limits^{12}_{j=1}H^{ji}_{12}q_j
	% H^{pq}_{ab}=1/2\sum^2_{j=1}L^p_{aj}L^p_{bj}
	% H^{pq}_{33}=1/2L^p_{33}L^q_{33}
	%
	% where 
	% L^i_{11} = dF_{1i}/dr
	% L^i_{22} = dF_{2i}/dz
	% L^i_{33} = F_{1i}/r
	% L^i_{12} = dF_{2i}/dr
	% L^i_{21} = dF_{1i}/dz
	% F_{ij} is a components of matrix of shape functions:
	%
	% F1 00	F2 00 F3 00 F4 00 F5 00 F6 00
	% 00 F1 00 F3 00 F3 00 F4 00 F5 00 F6
	
	L = l_matrix(el,r,z); % matrix of L components L11,L22,L33,L12,L21
	% L(1,i) means L^i_{11}
	% L(2,i) means L^i_{22}
	% L(3,i) means L^i_{33}
	% L(4,i) means L^i_{12}
	% L(5,i) means L^i_{21}

	for i = 1:12	
		for j = 1:12
			B(1,i) = B(1,i) + 0.5 * ( L(1,j)*L(1,i)+L(4,j)*L(4,i) )*q(j);
			B(2,i) = B(2,i) + 0.5 * ( L(5,j)*L(5,i)+L(2,j)*L(2,i) )*q(j);
			B(3,i) = B(3,i) + 0.5 * L(3,j)*L(3,i)*q(j);
			B(4,i) = B(4,i) + 0.5 * ( L(1,j)*L(5,i)+L(4,j)*L(2,i) )*q(j);
		end
	end


