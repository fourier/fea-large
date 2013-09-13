function [m,p] = apply_bc(M,P,boundary)
% function [m,p] = apply_bc(M,P,boundary)
% Apply fixed nodes boundary conditions to the stiffness matrix M and 
% right-side vector P
% Parameters:
% M - global stiffness matrix
% P - global vector of external loads
% boundary - matrix of boundary conditions

% constants
% EFree = 0 % free;
% EPrescribedX = 1 % : x prescribed 
% EPrescribedY = 2 % : y prescribed 
% EPrescribedXY = 3 % : x, y prescribed 
% EPrescribedZ = 4 % : z prescribed
% EPrescribedXZ = 5 % : x, z prescribed
% EPrescribedYZ = 6 % : y, z prescribed
% EPrescribedXYZ = 7  % : x, y, z prescribed.

	m = M;
	p = P;

	boundary_count = size(boundary,1);
	msize = length(m);

	% apply x-coordinate conditions
	for i = 1:boundary_count
		if contains_x(boundary(i,5))
			index = boundary(i,1)*3-2; % x-coordinate
			tmp = m(index,index);
			condition = boundary(i,2);
			for j = 1:msize
				p(j) = p(j) - m(index,j)*condition;
				m(index,j) = 0;
				m(j,index) = 0;
			end
			m(index,index) = tmp;
			p(index) = tmp*condition;
		end 
	end
	% apply y-coordinate conditions
	for i = 1:boundary_count
		if contains_y(boundary(i,5))
			index = boundary(i,1)*3-1; % x-coordinate
			tmp = m(index,index);
			condition = boundary(i,3);
			for j = 1:msize
				p(j) = p(j) - m(index,j)*condition;
				m(index,j) = 0;
				m(j,index) = 0;
			end
			m(index,index) = tmp;
			p(index) = tmp*condition;
		end 
	end
	% apply z-coordinate conditions
	for i = 1:boundary_count
		if contains_z(boundary(i,5))
			index = boundary(i,1)*3; % x-coordinate
			tmp = m(index,index);
			condition = boundary(i,4);
			for j = 1:msize
				p(j) = p(j) - m(index,j)*condition;
				m(index,j) = 0;
				m(j,index) = 0;
			end
			m(index,index) = tmp;
			p(index) = tmp*condition;
		end 
	end
end

function t = contains_x(idx)
  EFree = 0; % free;
  EPrescribedX = 1; % : x prescribed 
  EPrescribedY = 2; % : y prescribed 
  EPrescribedXY = 3; % : x, y prescribed 
  EPrescribedZ = 4; % : z prescribed
  EPrescribedXZ = 5; % : x, z prescribed
  EPrescribedYZ = 6; % : y, z prescribed
  EPrescribedXYZ = 7;  % : x, y, z prescribed.
	t = idx == EPrescribedX || idx == EPrescribedXY || ...
    idx == EPrescribedXZ || idx == EPrescribedXYZ ;
end    

function t = contains_y(idx) 
  EFree = 0; % free;
  EPrescribedX = 1; % : x prescribed 
  EPrescribedY = 2; % : y prescribed 
  EPrescribedXY = 3; % : x, y prescribed 
  EPrescribedZ = 4; % : z prescribed
  EPrescribedXZ = 5; % : x, z prescribed
  EPrescribedYZ = 6; % : y, z prescribed
  EPrescribedXYZ = 7;  % : x, y, z prescribed.
  t = idx == EPrescribedY || idx == EPrescribedXY || ...
    idx == EPrescribedYZ || idx == EPrescribedXYZ;
end

function t = contains_z(idx) 
  EFree = 0; % free;
  EPrescribedX = 1; % : x prescribed 
  EPrescribedY = 2; % : y prescribed 
  EPrescribedXY = 3; % : x, y prescribed 
  EPrescribedZ = 4; % : z prescribed
  EPrescribedXZ = 5; % : x, z prescribed
  EPrescribedYZ = 6; % : y, z prescribed
  EPrescribedXYZ = 7;  % : x, y, z prescribed.
  t = idx == EPrescribedZ || idx == EPrescribedXZ || ...
    idx == EPrescribedYZ || idx == EPrescribedXYZ;
end    

