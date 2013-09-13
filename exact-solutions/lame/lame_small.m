% Program for calculation Lame task for cylinder under
% small pressure

% Lame parameters
lambda = 100;
mu = 100;

% Geometry parameters
% height = 2*l
% inner radius = b
% outer radius = 2*a
b = 1; 
a = 2;
l = 8;


% Recalculate material parameters
G = mu;
nu = lambda/(2*(lambda+mu));
E = mu*(3*lambda+2*mu)/(lambda+mu);
% E = 2*G*(1+nu)

% Inner and outer pressure, Sigma_r(b) = -p0, Sigma_r(a) = -p1
p0 = 10000;
p1 = 0;

% Non-dimensional coordinates
x1 = b/a;
L = l/a;

% First case, not fixed z

% We want to calculate radial displacement at point r = 1,z = 1
r = 1;
z = 1;
x = r/a;
s = z/a;

Sigma_z0 = -2*nu*(p1*x1^2-p0)/(1-x1^2);
u0 = -(Sigma_z0*nu*a/E)*x;
w0 = (Sigma_z0*a/E)*s;

u = a/(2*G*(1-x1^2))*( ((1-nu)/(1+nu)) * (p1*x1^2-p0)*x + (p1-p0)*(x1^2)/x)
w = -(2*a*nu/(2*G*(1+nu))) *s* (p1*x1^2-p0)/(1-x1^2)

Sigma_r = (1/(1-x1^2)*x^2) * ((x^2-1)*p1*x1^2-(x^2-x1^2)*p0)
Sigma_z = (1/(1-x1^2)*x^2) * ((x^2+1)*p1*x1^2-(x^2+x1^2)*p0)