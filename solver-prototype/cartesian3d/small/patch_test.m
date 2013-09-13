function patch_test()
% patch test for elements 
% See Zienkiewitz v1, 6th edition, p.337-339
% using model A5 with coefficients: 
% l = 100, m = 100
%

% patch test using 1st element of the 'brick' geometry
% consider this element predefined
% How to get this element and therefore geometry
% el = element(nodes,elements,1)
% nodes = el
    
% Set up the element and nodes
    elements = [1,2,3,4,5,6,7,8,9,10];
    nodes = [0.5018, 2.4943, 0.5013;
             1.0000, 2.7482, 0.4974;
             1.0000, 2.5000, 1.0000;
             1.0000, 2.2394, 0.4979;
             0.7509, 2.6213, 0.4993;
             1.0000, 2.6241, 0.7487;
             0.7509, 2.4972, 0.7506;
             0.7509, 2.3668, 0.4996;
             1.0000, 2.4938, 0.4976;
             1.0000, 2.3697, 0.7489];
    
    % Set the solution considered
    % boundary - specify boundary conditions. it is a part of a
    % test
    X = test_C(nodes,elements);
    %X = test_C(nodes,elements)
    Y = exact_solution(nodes);
    X-exact_solution(nodes);

    %    [M,P] = apply_bc2(M,P,boundary);
    % Solve global system
    % X = inv(M)*P;
end

function X = test_A(nodes,elements)
    el = element(nodes,elements,1);
    K = local_stiffness(nodes);
    cond(K)
    u = exact_solution(nodes);
    X = K*u;
end

function X = test_B(nodes,elements)
    el = element(nodes,elements,1);
    K = local_stiffness(el);
    % add prescribed displacements to 5 nodes
    X = 0;
end

function X = test_C(nodes,elements)
    el = element(nodes,elements,1);
    K = local_stiffness(el);
    P = zeros(size(K,1),1);
    % prescribe displacements to 5 nodes: 
    % one node with all dofs and 4 nodes with 1 dof prescribed
    boundary = zeros(9,5);
    for i = 1:9
        boundary(i,1) = i;
    end
    % 9 fixed nodes:
    for i = 1:9
    boundary(i,5) = 7;
    boundary(i,2) = motion_x(nodes(i,1));
    boundary(i,3) = motion_y(nodes(i,2));
    boundary(i,4) = motion_z(nodes(i,3));
    end
    % nodes with prescribed displacements(let it be x-displacement)
    %  for i = 2:3
    %   boundary(i,5) = 1;
    %   boundary(i,2) = motion_x(nodes(i,1));
    %end
    [V,D] = eig(K);
    for i = 1:size(D)
        fprintf('%f\n',D(i,i))
    end
    [K,P] = apply_bc2(K,P,boundary);
    
    X = inv(K)*P;
end


function X = exact_solution(nodes)
    X = zeros(size(nodes,1)*3,1);
    for i = 1:size(nodes,1)
        X(3*i-2) = motion_x(nodes(i,1));
        X(3*i-1) = motion_y(nodes(i,2));
        X(3*i-0) = motion_z(nodes(i,3));
    end
end


function u = motion_x(x)
    u = 0.002*x;
end

function v = motion_y(y)
    v = -0.0006*y;
end

function w = motion_z(z)
    w = 0.000001*z;
end

function p = moved_node(nodes,i)
    p = [motion_x(nodes(i,1)),
         motion_y(nodes(i,2)),
         motion_z(nodes(i,3))];
end

% derivatives of the motion function
% s_ij = 1/2 * (du_i/du_j + du_j/du_i) 
function sx = strain_x(x,y,z)
    sx = 0.002;
end

function sy = strain_y(x,y,z)
    sy = -0.0006*y;
end

function sz = strain_z(x,y,z)
    sz = 0.000001;
end

function sxy = strain_xy(x,y,z)
    sxy = 0.;
end

function sxz = strain_xz(x,y,z)
    sxz = 0.;
end

function syz = strain_yz(x,y,z)
    syz = 0.;
end
