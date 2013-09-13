function h = h_function(e)
% Function for calculating h by given epsilon
% calculated using following maxima command:
%
% expand(determinant(matrix([2*e11+1,2*e12,0],[2*e12,2*e22+1,0],[0,0,2*e33+1])));
%

h = 4*e(1)*e(2)*e(3)+2*e(2)*e(3)-4*e(4)*e(4)*e(3)+2*e(1)*e(3)+...
		2*e(1)*e(2)-2*e(4)*e(4);	
