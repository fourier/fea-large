function result = qmr_analysis(flag,relres,iter,tab)
% Function qmr_analysis analyses results of QMR method and logs
% these results
% returns true if process was successfull or false otherwise
result = flag == 0;

switch flag
    case 0 
		string = sprintf('%sqmr converged to the desired tolerance %f within %d iterations',tab,relres,iter);
	case 1
		string = sprintf('%sqmr iterated maxit times but did not converge',tab,iter);
    case 2 
		string = sprintf('%sPreconditioner M was ill-conditioned',tab);
    case 3
    	string = sprintf('%sThe method stagnated. (Two consecutive iterates were the same)',tab);
    case 4
    	string = sprintf('%sOne of the scalar quantities calculated during qmr became too small or too large to continue computing.',tab);
	otherwise
		error('Wrong flag from QMR method received!');
end
logger(string);