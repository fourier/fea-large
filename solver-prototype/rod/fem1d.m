%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% MATLAB program for calculation 1d linear displacement
%% of solid rod under uniformed pressure causing displacement ux
%% Material model: Dimitrienko A5
%% Author: Veretennikov Alexey
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%% Constants definition
%%%%%%%%%%%%%%%%%%%%%%%
% geometry parameters
A = 1;
B = 7;
% material parameters
Lambda = 100;
Mu = 100;
E = Mu*(3*Lambda+2*Mu)/(Lambda+Mu);
% mesh parameters, number of elements
N = 3;
% boundary conditions
% displacement in last node
Ux = 0.05;
% structure of boundary conditions:
% x(0) = 0
% x(L) = Ux
%%%%%%%%%%%%%%%%%%
%% Mesh calculation
%%%%%%%%%%%%%%%%%%%
l = (B-A)/N;
nodes = (A:l:B)';
nodes_count = size(nodes,1);
elements = zeros(N,2);
for i = 1:N
	elements(i,1) = i;
	elements(i,2) = i+1;
end
elements_count = N;
element_size = size(elements,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparation of indexes etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize = nodes_count*1; % size of global matrix - number of nodes*dof
IDX = zeros(elements_count,2);
IDX = elements;
% fill boundary array;
boundary = zeros(2,2);
boundary(1,1) = 1; % index of node x(0)
boundary(1,2) = 0; % with prescribed displacement 0
boundary(2,1) = size(nodes,1); % index of nodex x(L)
boundary(2,2) = Ux; % with prescribed displacement Ux
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize global stiffness matrix and global force vector
M = zeros(msize,msize);
P = zeros(msize,1);
% loop over elements
for i = 1:elements_count
	% construct local stiffness matrix
	% sizes - number of dof per element
	local_size = element_size*1;
	K = ones(local_size,local_size); 
	K(1,2) = -1;
	K(2,1) = -1;
	S_i = nodes(elements(i,2))-nodes(elements(i,1));
	K = K*E/S_i;

	% distribute local in global matrix
	for k = 1:local_size
		for l = 1:local_size
			M(IDX(i,k),IDX(i,l)) = M(IDX(i,k),IDX(i,l)) + K(k,l);
		end
	end	
end	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(boundary,1)
	index = boundary(i,1);
	ux = boundary(i,2);
	tmp = M(index,index);
	for j = 1:msize
		P(j) = P(j) - M(index,j)*ux;
		M(index,j) = 0;
		M(j,index) = 0;
	end	
	M(index,index) = tmp;
	P(index) = tmp*ux;
end	
%%%%%%%%%%%%%%%%%%%%%%
%% Solve global system
%%%%%%%%%%%%%%%%%%%%%%
X = inv(M)*P;


