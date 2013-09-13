% Program for calculation uniaxial stresses from Neo-Hookean model

% Neo-Hookean model uses only 1 material parameter - shear modulus(G or mu)
global mu;
mu = 100;

% material model for Neo-Hookean material is T = -pE + mu*B
% where B is a Finger stress tensor, B = F*F^T
% taken from http://en.wikipedia.org/wiki/Neo-Hookean_solid

function t = Txx(k)
    global mu;
	t = mu*(k**2 - 1.0/k);
end

k = 1;
step = 0.008333;
for i=1:20
    stress = Txx(k);
    fprintf('%f %f\n',k,stress);
    k = k+step;
end
