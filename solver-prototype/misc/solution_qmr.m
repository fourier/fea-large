function [X,flag] = solution_qmr(A,b,tab)
% X - solution
% flag = true/false, if solution succeeded or not

%[L,U] = luinc(A,1e-14);
setup.type = 'crout';
setup.milu = 'row';
setup.droptol = 1e-14;
[L,U] = ilu(A,setup);
string = sprintf('%sCondition number = %f',tab,condest(A));
logger(string);
[X,flag,relres,iter] = qmr(A,b,1e-15,10000,L,U);
flag = qmr_analysis(flag,relres,iter,tab);

