function [nodes,elements,symmetry,boundary]=geometry
% Exports geometry of the task

%geometry_nodes;
%cylinder;
%cone;
%figure1;
dipl_K_4_ready;
% dipl_K_8_ready;

nodes = nds(:,2:3); 

elements = trngls+1;

symmetry = symm'+1;

boundary = bnd;
siz = size(boundary);
for i=1:siz(1)
	boundary(i,1) = boundary(i,1) + 1;
end
