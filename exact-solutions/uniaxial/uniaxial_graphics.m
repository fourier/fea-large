%
% Set constants
%

% starting strain
k1_begin = 1;
% step of strains
step = 0.008333;
% number of model
n = 5;

% material constants

%% For Yu.Dimitrienko task
%nu = 0.3;
%mu = 1;
%lambda = 2*mu*nu/(1-2*nu);

mu = 100;
lambda = 100;

betta = 100;
% number of iterations
iter = 100;

% prepare


% prepare file
fname = sprintf('uniaxial%d.txt',n);
f = fopen(fname,'w+');
% start loop
for i = 1:iter
	strain = k1_begin + i*step;
	[k2_an,s1_an,k2_bn,s2_bn] = uniaxial(strain,n,lambda,mu,betta);
%	string = sprintf("%f %f\n",strain,s1_an/(2*lambda*(1+nu)));
	string = sprintf('%f %f\n',(strain-1)*100,s1_an);
	fprintf(f,string);
end
fclose(f);
