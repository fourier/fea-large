function [nodes,elements,symmetry,boundary]=geometry_test
% Exports geometry of the task
nodes = [1, 1;
				 4, 1;
				 1, 4;
				 4, 4;
				 1, 7;
				 4, 7;
				 2.5, 1;
				 2.5, 2.5;
				 1, 2.5;
				 4, 2.5;
				 2.5, 4;
				 4, 5.5;
				 2.5, 7;
				 2.5, 5.5;
				 1, 5.5];

elements = [1,2,3,7,8,9;
						2,4,3,10,11,8;
						3,4,5,11,14,15;
						4,6,5,12,13,14];

boundary = [ 1, 0, 0; 
			  7, 0, 0;
			  2, 0, 0;
			  5, 0, 0.3;
				13,0, 0.3;
				6, 0, 0.3];	
symmetry = [];
