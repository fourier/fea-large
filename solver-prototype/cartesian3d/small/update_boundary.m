function [bnd] = update_boundary(nodes,boundary)
% [bnd] = update_boundary(nodes,boundary)
% Updates boundary conditions - make all nodes "free" in 
% all directions except Y and mark 1 point as fixed
    bnd = boundary;
    bnd(:,5) = 2;
    % find node with coordinates (1,1,1) and make it fixed	
    nodes_count = size(nodes,1);
    min_index = 1;
    max_index = 1;
    for i = 2:nodes_count
        if nodes(i,1) <= nodes(min_index,1) && ...
                nodes(i,2) <= nodes(min_index,2) && ...
                nodes(i,3) <= nodes(min_index,3)
            min_index = i;
        end
    end	
    for i = 2:nodes_count
        if nodes(i,1) == nodes(min_index,1) && ...
                nodes(i,2) >= nodes(max_index,2) && ...
                nodes(i,3) == nodes(min_index,3)
            max_index = i;
        end
    end	
    % set fixed boundary condition for this node
    for i = 1:size(boundary,1)
        if bnd(i,1) == min_index %|| ...
                                 %			bnd(i,1) == max_index
            bnd(i,5) = 7;
        end	
    end
    for i = 1:size(boundary,1)
        if bnd(i,3) == 0.05
             bnd(i,3) = -0.05;
        end
    end
    %    bnd(:,3) = 0.05;