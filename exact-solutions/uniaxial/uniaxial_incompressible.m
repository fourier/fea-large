% program for calculating incompressible uniaxial task for finite deformations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculation parameters and constants definitions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% Globals
%
% Lame constants
global lambda;
global mu;
global betta;
% Model type
global model_type;
% Model number
global model_number;

lambda = 100;
mu = 100;
betta = 100;

model_type = 'A';
model_number = 5;


function t = Txx(k)
    global model_type;
    if model_type == 'A' % model A(n)
        t = A_T2xx(k);
    else % model B(n)
        t = B_T2xx(k);
    end
end

function t = A_Txx(k)
    global model_number;
    global lambda;
    global mu;
    M = 0;
    N = 0;
    if model_number == 1
        M = 0.5*( 3 - 1/k**2 - 2*k )*(1.0/k**2 - k);
        N = (1.0/k**2) * ( 1 - 1./k**2 ) - k*(1-k);
    elseif model_number == 2
		M = 0.5*(3-1.0/k-2*sqrt(k))*(1.0/k-sqrt(k));
		N = 1.0/k*(1-1.0/k)-sqrt(k)*(1-sqrt(k));
    elseif model_number == 4
        M = 0.5*(k+2.0/sqrt(k)-3)*(k-1.0/sqrt(k));
        N = 1.0/sqrt(k)*(1-1.0/sqrt(k))-k*(1-k);
    elseif model_number == 5
        M = 0.5*(k**2+2.0/k-3)*(k**2-1.0/k);
        N = 1.0/k*(1-1.0/k)-k**2 *(1-k**2);
    end
    t = lambda*M+mu*N;
end

function t = A_T2xx(k)
    global model_number;
    global lambda;
    global mu;
    n = model_number;
    m = mu;
    l = lambda;
		kk1 = k^(n-3);
		kk2 = k^((3-n)/2.0);
    p = kk2*((2*l+2*m)*kk2+l*kk1-(3*l+2*m))/(n-3);
    t = -p + kk1*((l+2*m)*kk1+2*l*kk2-(3*l+2*m))/(n-3);
end


function t = B_Txx(k)
    global model_number;
    global mu;
    global betta;
    Z = 0;
    H = 0;
    if model_number == 1
        Z = 2*(k-1.0/(k**2));
    	H = k**2 - 1.0/k;
    elseif model_number == 2
		Z = sqrt(k)-1.0/k;
    	H = k - 1.0/sqrt(k);
    elseif model_number == 4
		Z = k - 1.0/sqrt(k);
		H = 1.0/k - sqrt(k);
    elseif model_number == 5
		Z = 2*(k**2-1.0/k);
		H = 1.0/(k**2) - k;
    end
    t = mu*(1+betta)*Z-mu*(1-betta)*H;			
end

% alternative implementation of Bn model stress
function t = B_T2xx(k) 
    global model_number;
    global mu;
    global betta;
    n = model_number;
    m = mu;
    b = betta;
    p = k^(0.5*(3-n))*m*(n-3)*(1+b+(1-b)*(k^(n-3)+k^(0.5*(3-n))));
    t = -p + (k^(n-3)) *m*(n-3)*(1+b+2*(1-b)*k^(0.5*(3-n)));
end

% k = (l+h)/l, where l is length of the truss
% 1 means what truss is fixed
k = 1;
step = 0.008333;

for i=1:20
    stress = Txx(k);
    fprintf('%f %f\n',k,stress);
    k = k+step;
end
