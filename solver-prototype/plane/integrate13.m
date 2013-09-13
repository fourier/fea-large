function M = integrate13(el,F)
% Function for integrating matrixes on element
% by given element el, pointer to matrix function @F 
% and vector of solutions or any other argument q
% Integration formulae by 13 points from lecture from 
% University of Tokyo, Nonlinear Finit Element Analysis
% Lecture 3
M = [];
rows = 0;
cols = 0;

siz = size(F(el,1,1));
rows = siz(1);
cols = siz(2);

% calculate nodes for integration
% r,s - coordinates, w - weights
r = zeros(13,1);
s = zeros(13,1);
w = zeros(13,1);

r(1)  = 0.0651301029002; s(1)  = r(1); 							w(1)  = 0.0533472356008;
r(2)  = 0.8697397941956; s(2)  = r(1); 							w(2)  = w(1);
r(3) 	= r(1);						 s(3)  = r(2); 							w(3)  = w(1);
r(4) 	= 0.3128654960049; s(4)  = 0.0486903154253;		w(4)  = 0.0771137608903;
r(5) 	= 0.6384441885698; s(5)  = r(4);							w(5)  = w(4);
r(6) 	= s(4);						 s(6)  = r(5);							w(6)  = w(4);
r(7) 	= r(5);						 s(7)  = r(6);							w(7)  = w(4);
r(8) 	= r(4);						 s(8)  = r(5);							w(8)  = w(4);
r(9) 	= r(6);						 s(9)  = r(4);							w(9)  = w(4);
r(10) = 0.2603459660790; s(10) = r(10);							w(10) = 0.1756152574332;
r(11)	= 0.4793080678419; s(11) = r(10);							w(11) = w(10);
r(12) = r(10);					 s(12) = r(11);							w(12) = w(10);
r(13) = 1/3.0;					 s(13) = r(13);							w(13) = -0.149570044467;

p = zeros(13,2);

for i=1:13
	p(i,1) = el(1,1)*r(i)+el(2,1)*s(i)+el(3,1)*(1-r(i)-s(i));
	p(i,2) = el(1,2)*r(i)+el(2,2)*s(i)+el(3,2)*(1-r(i)-s(i));
end

% ok, let's do the integration by using formulae
% \int(F)dS = \sum_{i=1}^{13} f(x_i)*w_i

M = zeros(rows,cols);
for i = 1:13
	M = M + w(i)*F(el,p(i,1),p(i,2))*p(i,1);		
end	
