% Program for calculation uniaxial stresses from Chernyh model

% Material properties of Chernyh model
global mu;
global betta;
mu = 100;
betta = 100;

function t = Txx(k)
    global mu;
	global betta;
	t = mu*((1+betta)*(k-1.0/sqrt(k))+(1-betta)*(sqrt(k)-1.0/k));
end

k = 1;
step = 0.008333;
for i=1:100
    stress = Txx(k);
    fprintf('%f %f\n',k,stress);
    k = k+step;
end
