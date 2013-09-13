function [a,v] = create_bc(A,V,idx,val)
    [a,v] = create_bc2(A,V,idx,val);
end

function [a,v] = create_bc1(A,V,idx,val)
% works properly
    a = A;
    v = V;
    msize = size(A,1);
    tmp = a(idx,idx);
    for j = 1:msize
        v(j) = v(j) - a(j,idx)*val;
        a(idx,j) = 0;
        a(j,idx) = 0;
    end
    a(idx,idx) = tmp;
    v(idx) = tmp*val;
end

function [a,v] = create_bc2(A,V,idx,val)
% works properly
    a = A;
    v = V;
    msize = size(A,1);
    for j = 1:msize;
        a(idx,j) = 0;
    end
    a(idx,idx) = 1;
    v(idx) = val;
end

function [a,v] = create_bc3(A,V,idx,val)
% works properly
    a = A;
    v = V;
    Alpha = 1e8;
    new_k = a(idx,idx)*Alpha;
    a(idx,idx) = new_k;
    v(idx) = val*new_k;
end


