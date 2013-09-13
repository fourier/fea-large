% program for calculating incompressible Lame task for finite deformations

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
% Geometry constants
% inner radius
global R1;
% outer radius 
global R2;
% Model type
global model_type;
% Model number
global model_number;
% Precision
global eps;
% Inner pressure
global Q1;
% Outer pressure
global Q2;
% k = (l+h)/l, where l is length of the cylinder
% 1 means what cylinder is fixed
global K;

% Set global values
lambda = 100;
mu = 1;
betta = 1;
R1 = 1;
R2 = 2;
model_type = 'B';
model_number = 5;
eps = 1e-5;
Q1 = 0;
Q2 = 0;
K = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Auxulary functions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function for calculation of the integral of function of 1st order
% with the name stored in variable func_name with arguments param1 
% and param2 with limits a and b
% function call signature:
% func_name(x,param1,param2)
function result = simps_integrate(func_name,param1,param2,a,b)
    % number of points in numeric integration + 1
    count = 101;
    n = (count - 1)/2;
    result = 0;
    H = (b-a)/(6.*n);
    lim = 2*n;
    for i=1:lim-1
        x = (b-a)*i/(count-1) + a;
        if mod(i,2) == 1
            result = result + 4*feval(func_name,x,param1,param2);
        else
            result = result + 2*feval(func_name,x,param1,param2);
        end
    end
    result = result + feval(func_name,a,param1,param2) + ...
        feval(func_name,b,param1,param2);		
    result = result*H;
end

function r = test_func(x,arg1,arg2)
    r = x*x*sin(x)+arg1*arg2;
end

% test integration:
% a = simps_integrate("test_func",4,6,2,10)

% blowing function, or new radius
% in point r
function r = f(r,C,k)
    r = sqrt(r*r/k + C);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Elastic Stresses
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1st invariant of the right Cauchy deformation tensor C(n)
function i = I1_C(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    i = ( ( r/(F*k) ) ** (n-3) + ...
          ( F/r ) ** (n-3) + k ** (n-3) - 3 )/(n-3);
end

%%
%% Cauchy stresses
%%

%%
%% Model A
%%

% r component of elastic part of Cauchy stress tensor for model A(n)
function t = A_Ts_rr(r,C,k)
    global model_number;
    global lambda;
    global mu;
    n = model_number;
    F = f(r,C,k);
    t = lambda*I1_C(r,C,k)+2*mu*( (r/(F*k)) ** (n-3) - 1 )/(n-3);
end

% phi component of elastic part of Cauchy stress tensor for model A(n)
function t = A_Ts_ff(r,C,k)
    global model_number;
    global lambda;
    global mu;
    n = model_number;
    F = f(r,C,k);
    t = lambda*I1_C(r,C,k)+2*mu*( (F/r) ** (n-3) - 1 )/(n-3);

end

% z component of elastic part of Cauchy stress tensor for model A(n)
function t = A_Ts_zz(r,C,k)
    global model_number;
    global lambda;
    global mu;
    n = model_number;
    F = f(r,C,k);
    t = lambda*I1_C(r,C,k)+2*mu*( k ** (n-3) - 1 )/(n-3);
end


%%
%% Model B
%%

% r component of elastic part of Cauchy stress tensor for model B(n)
function t = B_Ts_rr(r,C,k)
    global model_number;
    global betta;
    global mu;
    n = model_number;
    F = f(r,C,k);
    m = mu;
    b = betta;
    t = m*(n-3)* ( 1+b + (1-b)*((F/r) ** (n-3) + k**(n-3)) );
end

% phi component of elastic part of Cauchy stress tensor for model B(n)
function t = B_Ts_ff(r,C,k)
    global model_number;
    global betta;
    global mu;
    n = model_number;
    F = f(r,C,k);
    m = mu;
    b = betta;
    t = m*(n-3)* ( 1+b + (1-b)*( (r/(F*k)) ** (n-3) + k**(n-3)) );
end


% z component of elastic part of Cauchy stress tensor for model B(n)
function t = B_Ts_zz(r,C,k)
    global model_number;
    global betta;
    global mu;
    n = model_number;
    F = f(r,C,k);
    m = mu;
    b = betta;
    t = m*(n-3)* ( 1+b + (1-b)*( (F/r) ** (n-3) + (r/(F*k)) ** (n-3) ) );
end


%%
%% Model independent stress
%%

% r-component of stress tensor T
function t = Ts_rr(r,C,k)
    global model_type;
    if model_type == 'A' % model A(n)
        t = A_Ts_rr(r,C,k);
    else % model B(n)
        t = B_Ts_rr(r,C,k);
    end
end

% phi-component of stress tensor T
function t = Ts_ff(r,C,k)
    global model_type;
    if model_type == 'A' % model A(n)
        t = A_Ts_ff(r,C,k);
    else % model B(n)
        t = B_Ts_ff(r,C,k);
    end
end

% z-component of stress tensor T
function t = Ts_zz(r,C,k)
    global model_type;
    if model_type == 'A' % model A(n)
        t = A_Ts_zz(r,C,k);
    else % model B(n)
        t = B_Ts_zz(r,C,k);
    end
end


%%
%% Components of true Cauchy stress tensor
%%

% r component of true Cauchy stress tensor
function t = T_rr(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = (r/(F*k))**(n-3); 
    t = temp*Ts_rr(r,C,k);
    t = -hydrostatic_pressure(r,C,k)+t;
end

% phi component of true Cauchy stress tensor
function t = T_ff(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = (F/r)**(n-3); 
    t = temp*Ts_ff(r,C,k);
    t = -hydrostatic_pressure(r,C,k)+t;
end

% z component of true Cauchy stress tensor
function t = T_zz(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = k**(n-3); 
    t = temp*Ts_zz(r,C,k);
    t = -hydrostatic_pressure(r,C,k)+t;
end


%%
%% 1st Piola stresses
%%

% r component of elastic part of Piola stress tensor 
function p = Ps_rr(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = (r/(F*k))**(n-3-1); 
    p = temp * Ts_rr(r,C,k);
end

% phi component of elastic part of Piola stress tensor 
function p = Ps_ff(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = (F/r)**(n-3-1); 
    p = temp * Ts_ff(r,C,k);
end

% z component of elastic part of Piola stress tensor 
function p = Ps_zz(r,C,k)
    global model_number;
    n = model_number;
    F = f(r,C,k);
    temp  = k**(n-3-1); 
    p = temp * Ts_zz(r,C,k);
end

%%
%% h-function
%%

% Derivative of Piola stress compnent rr by r coordinate
function d = dPs_rr_dr(r,C,k)
    global R1;
    global R2;
    dr = (R2-R1)/1000.;
    d = (Ps_rr(r+dr,C,k)-Ps_rr(r,C,k))/dr;
end

function H = h_dimitrienko(r,C,k)
  global model_number;
  n = model_number;
  global R1;
  global R2;
  dr = (R2-R1)/1000.;
  f1 = f(r,C,k);
  f2 = f(r+dr,C,k);
  Sigma1 = f1*k*Ts_rr(r,C,k)*(r/(f1*k))^(n-3);
  Sigma2 = f2*k*Ts_rr(r+dr,C,k)*((r+dr)/(f2*k))^(n-3);
  dSigma = (Sigma2-Sigma1)/dr;
  Sigma3 = Ts_ff(r,C,k)*(f1/r)^(n-3);
  H = dSigma - Sigma3*r/f1;
end

% h function itself
function H = h(r,C,k)
  H = r*dPs_rr_dr(r,C,k)+Ps_rr(r,C,k)-Ps_ff(r,C,k);
%  H = h_dimitrienko(r,C,k);
end

% function under integral
function u = under_integral(r,C,k)
    u = h(r,C,k)/f(r,C,k);
end

%%
%% Hydrostatic pressure
%% 

% P0 constant
function P0 = p0(r,C,k)
    % import globals
    global R1;
    global Q1;
    % auxulary calculations
    f1 = f(R1,C,k);
    P1 = Ps_rr(R1,C,k);
    P0 = (P1 + Q1)*R1/(f1*k);
end

% p(r) - hydrostatic pressure
function P = hydrostatic_pressure(r,C,k)
    P = 0;
    % import globals
    global R1;
    % auxulary calculations
    integr = simps_integrate("under_integral",C,k,R1,r);
    % calculation of a function  
    P = p0(r,C,k) + 1/k*integr;
end

%%
%% Function to find roots for
%% F(C) = 0
%% with boundary conditions in reference configuration
%%
function F = functional_C_piola(C)  
    % import globals
    global R1;
    global R2;
    global K;
    global Q1;
    global Q2;
    % auxulary calculations
    f1 = f(R1,C,K);
    f2 = f(R2,C,K);
    P1 = Ps_rr(R1,C,K);
    P2 = Ps_rr(R2,C,K);
    integr = simps_integrate("under_integral",C,K,R1,R2);
    % calculation of a function  
    F = -R1/R2*(f2/f1)*(P1+Q1) + 1/K*integr+P2+Q2;
end

%%
%% Function to find roots for
%% F(C) = 0
%% with boundary conditions in actual configuration
%%
function F = functional_C(C)  
    % import globals
    global model_number;
    global R1;
    global R2;
    global K;
    global Q1;
    global Q2;
    % auxulary calculations
    n = model_number;
    f1 = f(R1,C,K);
    f2 = f(R2,C,K);
    Sigma1 = Ts_rr(R1,C,K)*(R1/(f1*K))^(n-3);
    Sigma2 = Ts_rr(R2,C,K)*(R2/(f2*K))^(n-3);
    integr = simps_integrate("under_integral",C,K,R1,R2);
    % calculation of a function  
    F = Sigma1-Sigma2+(Q1-Q2)+1/K*integr;
end


function calcc()
  global model_type;
  global model_number;
  global R1;
  global R2;
  global K;
  global Q1;
  global Q2;
  global mu;
  global betta;


  model_type = 'B';
  mu = 1;
  models = [1,2,4,5];
  betas = [-0.9, -0.5, 0, 0.5, 1];
  %Q1s = [0.00001, 0.00005, 0.0001];
  Q1s = [0.001, 0.005, 0.01];

  Q2 = 0;
  Q1 = 0;
  for q1 = Q1s
    Q1 = q1;
    for b = betas
      betta = b;
      for i = models
        model_number = i;
        [C,info] = fsolve(@functional_C,0);
        fname = sprintf('plots/model%c%d_beta%.2f_Q1is%.5f_Q2is%.2f_Srr.txt',model_type,model_number,betta,Q1,Q2);
        f = fopen(fname,'w+');
        steps = 20;
        h = (R2-R1)/steps; 
        for i = 0:steps
          point = R1+h*i;
          stress1 = T_ff(point,C,K);
          fprintf(f,'%e %e\n',point,stress1/mu);
        end
        fclose(f);
      end
    end
  end
end


function plot_c()
  model_type = 'B';
  mu = 1;
  models = [1,2,4,5];
  betas = [-0.9, -0.5, 0, 0.5, 1];
  Q1s = [0.00001, 0.0005, 0.0001];

  K = 1;
  Q2 = 0;
  Q1 = 0;

  for q1 = Q1s
    Q1 = q1;
    for b = betas
      betta = b;
      for i = models
        model_number = i;
        fname = sprintf('plots/q1_is_%f_beta_is_%f_n%d',Q1,betta,model_number);
        f = fopen(fname,'w+');
        a = -1.0;
        b = 1.;
        steps = 100;
        for i = 1:steps
          point = a+i*(b-a)/(steps+1);
          fprintf(f,'%e %e\n',point,functional_C(point));
        end 
        fclose(f);
      end
    end
  end  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
%% Calculation 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc(steps)
  models = [1,2,4,5];
  betas = [-0.5, 0, 0.5, 1];

  for b = betas
    betta = b;
    for i = models
      model_number = i;
      fname = sprintf('plots/lame_dQ_Srr_beta%.1f_model_%s%d.txt',betta,model_type,model_number);
      f0 = fopen(fname,'w+');

      fname1 = sprintf('plots/lame_dQ_Fc_beta%.1f_model_%s%d.txt',betta,model_type,model_number);
      f1 = fopen(fname1,'w+');
      fname2 = sprintf('plots/lame_dQ_radial_def_beta%.1f_model_%s%d.txt',betta,model_type,model_number);
      f2 = fopen(fname2,'w+');

      K = 1;
      Q1 = 0;
      step = -0.05;
      for i=1:steps
        Q1 = Q1 + step;
          [C,info] = fsolve(@functional_C,0);
          point = (R2+R1)/2.0;
          stress = T_rr(point,C,K)/mu;
        radius_int = f(R1,C,K);
        radius_ext = f(R2,C,K);
        thickness = radius_ext-radius_int;
        old_thickness = R2-R1;
        deformation = 100*(thickness-old_thickness)/old_thickness;
          fprintf(f0,'%f %f\n',Q1,stress);
        fprintf(f1,'%f %f\n',Q1,C);
        fprintf(f2,'%f %f\n',Q1,deformation);
      end
      fclose(f0);
      fclose(f1);
      fclose(f2);
    end
  end
end

function check_with_uniaxial(steps)
  global mu;
  global betta;
  global Q1;
  global Q2;
  global K;
  global model_type;
  global model_number;
  global R1;
  global R2;

  mu = 1;
  betta = 1;
  Q1 = 0;
  Q2 = 0;
  K = 1;
  model_type = 'B';
  model_number = 5;

  step = 0.008333;
  for i=1:20
      [C,info] = fsolve(@functional_C,0);
    point = (R2+R1)/2.0;
      stress = T_zz(point,C,K);
      fprintf('%f %f\n',K,stress);
      K = K+step;
  end

end

function dq_ds()

  global mu;
  global betta;
  global Q1;
  global Q2;
  global K;
  global model_type;
  global model_number;
  global R1;
  global R2;

model_type = 'B';
betta = 1;
R1 = 1;
R2 = 1.01;

models = [1];%,2,4,5];
for i = models
	model_number = i;
		fname = sprintf('plots/dQ_dS_model_%s%d.txt',model_type,model_number);
		f0 = fopen(fname,'w+');
		K = 1;
		Q1 = 0;
		stepp = 0.1;
    steps = 50;
		for i=1:steps
			Q1 = Q1 + stepp;
		  [C,info] = fsolve(@functional_C,0);
%		  point = (R2+R1)/2.0;
			radius_int = f(R1,C,K);
			radius_ext = f(R2,C,K);
%			thickness = radius_ext-radius_int;
%			old_thickness = R2-R1;
%			deformation = 100*(thickness-old_thickness)/old_thickness;
      dr = (R2-R1)/100.0;
      deformation = 100+(radius_int-R1)/dr;
			fprintf(f0,'%f %f\n',Q1,deformation);
%      fprintf(f0,'%f %f\n',radius_int,Q1);
%      fprintf(f0,'%f %e\n',Q1,C);
		end
		fclose(f0);
end

end
%check_with_uniaxial(20);
%calc(20)
calcc()