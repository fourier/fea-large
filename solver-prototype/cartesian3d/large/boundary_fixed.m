function [bnd] = boundary_fixed(boundary)
% function [bnd] = boundary_fixed(boundary)
% Updates boundary conditions - only "fixed" nodes remain

    count_fixed = 0;
    for i = 1:size(boundary,1)
        if boundary(i,2) == boundary(i,3) == boundary(i,4) == 0
            count_fixed = count_fixed + 1;
        end
    end
    bnd = zeros(count_fixed,5);
    index = 1;
    for i = 1:size(boundary,1)
        if boundary(i,2) == boundary(i,3) == boundary(i,4) == 0
            bnd(index,:) = boundary(i,:);
            index = index + 1;
        end
    end
    
