function [nodes,elements,symmetry,boundary]=geometry_test
% Exports geometry of the task
nodes = [2,2;
				 5,2;
				 2,5;
				 5,5;
				 3.5,2;
				 3.5,3.5;
				 2,3.5;
				 5,3.5;
				 3.5,5];


elements = [1,2,3,5,6,7;
						2,4,3,8,9,6];

boundary = [ 1,0,0;
						 5,0,0;
						 2,0,0;
						 3,0,0.1;
						 9,0,0.1;
						 4,0,0.1];
symmetry = [];
