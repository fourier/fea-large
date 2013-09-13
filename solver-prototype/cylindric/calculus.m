% Calc task of incompressible/cvompressible finite deformations of 
% cylindric bodies

addpath ('6p3_cylindric','elements','geometries');
addpath ('linear','nonlinear','condder');
addpath ('misc','methods');
addpath ('materials');
addpath ('tests');

% Name of file for logging
global g_logfilename;
% Handle to open file for logging
global g_logfile;
% Defines if we will calculate a compressible of incompressible finite deformations
global g_compressible;
% Path where we will store results of calculations
global g_resultspath;
% Should we use Jacobi refinement in solution
global g_jacobi_refinement;
% Should we use line search refinement in solution
global g_linesearch;
% Model: 'A5' or 'neohookean-compressible'
global g_model;
% material constants
global g_lambda;
global g_mu;


g_logfilename = 'results/calculus.log';
g_logfile = fopen(g_logfilename,'w+');
g_resultspath = 'results';
g_compressible = true;
g_jacobi_refinement = false;
g_linesearch = false;
g_model = 'A5';
%g_model = 'neohookean-compressible';
%g_lambda = 2.74;
%g_mu = 0.69;
g_lambda = 100;
g_mu = 100;

mkdir 'results';

% Start calculations
%feaincrement;
%fealinear;
conditionderiv;

% Close log file
fclose(g_logfile);
clear g_logfile;
