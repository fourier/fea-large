function M = integrate7(el,F)
% Function for integrating matrixes on element
% by given element el, pointer to matrix function @F 
% and vector of solutions q
% Integration formulae by 7 points from Zienkiewicz, vol1, p 222
M = [];
rows = 0;
cols = 0;
siz = size(F(el,1,1));
rows = siz(1);
cols = siz(2);

% coefficients
A1 = 0.0597158717;
B1 = 0.4701420641;
A2 = 0.7974269853;
B2 = 0.1012865073;

% weights 
w1 = 0.225;
w2 = 0.1323941527;
w3 = 0.1259391805;

% calculate nodes for integration
elem = el(1:3,:);

a = [ 1/3.0, 1/3.0, 1/3.0 ]*elem;
b = [ A1, B1, B1 ]*elem;
c = [ B1, A1, B1 ]*elem;
d = [ B1, B1, A1 ]*elem;
e = [ A2, B2, B2 ]*elem;
f = [ B2, A2, B2 ]*elem;
g = [ B2, B2, A2 ]*elem;

% ok, let's do the integration by using formulae
% \int(F)dS = \sum_{i=1}^7 f(x_i)*w_i

M = w1*F(el,a(1),a(2))*a(1)+...
		w2*F(el,b(1),b(2))*b(1)+w2*F(el,c(1),c(2))*c(1)+...
		w2*F(el,d(1),d(2))*d(1)+...
		w3*F(el,e(1),e(2))*e(1)+w3*F(el,f(1),f(2))*f(1)+w3*F(el,g(1),g(2))*g(1);

