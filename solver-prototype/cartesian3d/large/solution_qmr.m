function [X,flag] = solution_qmr(A,b)
% X - solution
% flag = true/false, if solution succeeded or not

[L,U] = luinc(A,1e-14);
[X,flag,relres,iter] = qmr(A,b,1e-15,10000,L,U);
flag = qmr_analysis(flag,relres,iter);
