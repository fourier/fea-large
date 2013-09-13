% Test for the external pressure forces

nodes = [
    0,   0,   0;
    1,   0,   0;
    0,   1,   0;
    0,   0,   1;
    1/2, 0,   0;
    1/2, 1/2, 0;
    0,   1/2, 0;
    0,   0,   1/2;
    1/2, 0,   1/2;
    0,   1/2, 1/2;
];

elements = [
    1,2,3,4,5,6,7,8,9,10;
];

el = element(nodes,elements,1);

% Sides of the triangle, node numbers:
side1 = [1,2,3,5,6,7];
side2 = [2,3,4,6,10,9];
side3 = [1,2,4,5,9,8];
side4 = [1,3,4,7,10,8];

dxis = zeros(3,2);

for i = 1:3
    for j = 1:2
        for k = 1:6
        end 
    end
end