function [x,i] = cg(A,b,x0)
  x = x0;
  r = b - A*x0;
  p = r;
  for i = 1:10000
    temp = A*p;
    a1 = r'*r;
    a2 = (temp'*p);
    alpha = a1/a2;
    x = x + alpha*p;
    r = r - alpha*temp;
    if max(abs(r)) < 1e-15
      break;
    end
    a2 = r'*r;
    beta = a2/a1;           
    p = r + beta*p;
  end
end