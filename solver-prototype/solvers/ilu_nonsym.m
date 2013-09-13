function [L,U] = ilu(A)
  L = zeros(size(A));
  U = zeros(size(A));
  n = size(A,1);
  
  for k = 1:n
    for j = 1:k-1
      if A(k,j) ~= 0
        sum = 0;
        for i = 1:j-1
          sum = sum + L(k,i)*U(i,j);
        end
        L(k,j) = (A(k,j) - sum)/U(j,j);
      end
    end
    L(k,k) = 1;
    for j = k:n
      if A(k,j) ~= 0
        sum = 0;
        for i = 1:k-1
          sum = sum + L(k,i)*U(i,j);
        end
        U(k,j) = A(k,j)-sum;
      end
    end
  end
end