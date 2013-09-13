function [x,i] = qmr(A,b,x0)
  x = x0;
  r = b - A*x0;
  gamma = norm(r);
  w = v = r/gamma;
  for i = 1:1
    [v,w,alpha,delta] = lanzos_biorthogonalization(A,v,w);
    

    
  end
end


function [v,w,a,d] = lanzos_biorthogonalization(A,v,w)
  b = d = 0;
  w0 = v0 = zeros(size(v));
  for i = 1:size(A,1)
    a = (A*v)'*w;
    v1 = A*v-a*v-b*v0;
    w1 = A'*w-a*w-d*w0;
    d = sqrt(v1'*w1);
    if d < 1e-15
      break;
    end
    b = v1'*w1/d;
    w0 = w;
    w = w1/b;
    v0 = v;
    v = v1/d;
  end
  w'*v
end  