function [X,flag] = solution_pcg(A,b,tab)
% X - solution
% flag = true/false, if solution succeeded or not

[L,U] = luinc(A,1e-14);
string = sprintf('%sCondition number = %f',tab,condest(A));
logger(string);
[X,flag,relres,iter] = pcg(A,b,1e-15,10000,L,U);
flag = pcg_analysis(flag,relres,iter,tab);

