function M = integrate7initial(el0,el,F)
% function M = integrate7initial(el0,el,F)
% Function for integrating matrixes on element
% by given elements el0,el and pointer to matrix function @F 
% Integration formulae by 7 points from Zienkiewicz, vol1, p 222
M = [];
cntr = center_elem(el0);
siz = size(F(el0,el,cntr(1),cntr(2),cntr(1),cntr(2)));
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
elem0 = el0(1:3,:);
elem = el(1:3,:);

% initial configuration nodes
a0 = [ 1/3.0, 1/3.0, 1/3.0 ]*elem0;
b0 = [ A1, B1, B1 ]*elem0;
c0 = [ B1, A1, B1 ]*elem0;
d0 = [ B1, B1, A1 ]*elem0;
e0 = [ A2, B2, B2 ]*elem0;
f0 = [ B2, A2, B2 ]*elem0;
g0 = [ B2, B2, A2 ]*elem0;
% actual configuration nodes
a = [ 1/3.0, 1/3.0, 1/3.0 ]*elem;
b = [ A1, B1, B1 ]*elem;
c = [ B1, A1, B1 ]*elem;
d = [ B1, B1, A1 ]*elem;
e = [ A2, B2, B2 ]*elem;
f = [ B2, A2, B2 ]*elem;
g = [ B2, B2, A2 ]*elem;



% ok, let's do the integration by using formulae
% \int(F)dS = \sum_{i=1}^7 f(x_i)*w_i
M = w1*F(el0,el,a0(1),a0(2),a(1),a(2))*a(1)+...
		w2*F(el0,el,b0(1),b0(2),b(1),b(2))*b(1)+...
		w2*F(el0,el,c0(1),c0(2),c(1),c(2))*c(1)+...
		w2*F(el0,el,d0(1),d0(2),d(1),d(2))*d(1)+...
		w3*F(el0,el,e0(1),e0(2),e(1),e(2))*e(1)+...
		w3*F(el0,el,f0(1),f0(2),f(1),f(2))*f(1)+...
		w3*F(el0,el,g0(1),g0(2),g(1),g(2))*g(1);
