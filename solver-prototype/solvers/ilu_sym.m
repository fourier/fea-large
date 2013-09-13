function [L,D] = ilu_sym(A)
  L = zeros(size(A));
  D = zeros(size(A,1),1);
  n = size(A,1);
  
  for k = 1:n
    for j = 1:k-1
      if A(k,j) ~= 0
        sum = 0;
        for i = 1:j-1
          sum = sum + D(i)*L(i,k)*L(j,k);
        end
        L(k,j) = (A(k,j) - sum)/D(j);
      end
    end
    sum = 0;
    L(k,k) = 1;
    for i = 1:k-1
      sum = sum + (L(k,i)^2)*(D(i)^2);
    end
    D(k) = A(k,k)-sum;
  end
end