function N = construct_n(el,q,r,z)
% Function for constructing N matrix
% by given element el and point (r,z)

B = b_matrix(el,r,z);
B1 = b1_matrix(el,q,r,z);
B2 = b2_matrix(el,q,r,z);
Epsilon = epsilon_nonlin(el,q,r,z);
R = r_matrix(Epsilon);
Fp = fp_matrix(el,r,z);

I0 = [1;1;1;0];
I = ones(4);
I(4,:) = zeros(1,4);
I(:,4) = zeros(4,1);
I = I - eye(4); 

N = (B'+B1'+B2')*(I0+2*(I+R)*(B+B1)*q)*Fp';

