function result = pcg_analysis(flag,relres,iter,tab)
% Function qmr_analysis analyses results of QMR method and logs
% these results
% returns true if process was successfull or false otherwise
result = flag == 0;

switch flag
    case 0 
		string = sprintf('%spcg converged to the desired tolerance %f within %d iterations',tab,relres,iter);
	case 1
		string = sprintf('%sqmr iterated maxit times but did not converge',tab,iter);
    case 2 
		string = sprintf('%sPreconditioner M was ill-conditioned',tab);
    case 3
    	string = sprintf('%sPCG method stagnated. (Two consecutive iterates were the same)',tab);
    case 4
    	string = sprintf('%sOne of the scalar quantities calculated during PCG became too small or too large to continue computing.',tab);
	otherwise
		error('Wrong flag from PCG method received!');
end
logger(string);