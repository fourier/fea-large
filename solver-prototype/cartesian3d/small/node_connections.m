function conn = node_connections(idx,nodes,elements)

  %% find nodal connections for the node with index idx
  arr = zeros(1,size(nodes,1));
  index = 0;
  for i = 1:size(elements,1) % loop by elements
    for j = 1:size(elements,2) % loop by nodes per element
      if elements(i,j) == idx
        if j == 1 % our node is the first in element
          index = index + 1;
          arr(index) = elements(i,j+1);
        elseif j == size(elements,2) % our node is the last in element
          index = index + 1;
          arr(index) = elements(i,j-1);
        else % otherwise our node somewhere in the center of element
          index = index + 1;
          arr(index) = elements(i,j-1);
          index = index + 1;
          arr(index) = elements(i,j+1);
        end
      end
    end
  end
  
  arr = unique(arr);
  conn = arr(2:size(arr,2));

end