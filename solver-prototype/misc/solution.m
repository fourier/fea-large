function [X,flag] = solution(A,b,tab)
% X - solution
% flag = true/false, if solution succeeded or not
 X = inv(A)*b;
 flag = true;
%[X,flag] = solution_qmr(A,b,tab);
%[X,flag] = solution_pcg(A,b,tab);
