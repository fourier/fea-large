function result = qmr_analysis(flag,relres,iter)
% Function qmr_analysis analyses results of QMR method and logs
% these results
% returns true if process was successfull or false otherwise
result = flag == 0;

switch flag
    case 0 
		string = sprintf('qmr converged to the desired tolerance %f within %d iterations',relres,iter);
	case 1
		string = sprintf('qmr iterated maxit times but did not converge',iter);
    case 2 
		string = sprintf('Preconditioner M was ill-conditioned');
    case 3
    	string = sprintf('The method stagnated. (Two consecutive iterates were the same)');
    case 4
    	string = sprintf('One of the scalar quantities calculated during qmr became too small or too large to continue computing.');
	otherwise
		error('Wrong flag from QMR method received!');
end
fprintf(string);