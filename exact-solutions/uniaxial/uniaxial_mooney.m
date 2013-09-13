% Program for calculation uniaxial stresses from Mooney-Rivelin model

% Material properties of Mooney-Rivelin model
global mu;
global lambda;
mu = 100;
lambda = 100;

% material model for Neo-Hookean material is T = -pE + 2*C1*B+2*C2*B^-1
% where B is a Finger stress tensor, B = F*F^T
% C1 = 0.5*G or 0.5*mu, where G (or mu) - shear modulus
% taken from http://en.wikipedia.org/wiki/Mooney-Rivlin_solid

function t = Txx(k)
    global mu;
	global lambda;
	C1 = 0.5*mu;
	C2 = lambda;
	t = (2*C1-2*C2/k)*(k**2-1.0/k);
end

k = 1;
step = 0.008333;
for i=1:20
    stress = Txx(k);
    fprintf('%f %f\n',k,stress);
    k = k+step;
end
