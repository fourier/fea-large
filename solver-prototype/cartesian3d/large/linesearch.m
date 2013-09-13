function eta = linesearch(R,X,eta,r0,r,desired_tolerance)
% function eta = linesearch(R,X,eta,r0,r)
    r1 = X'*R;
    if abs(r1) < desired_tolerance
        return
    end
    alpha = r0/r
    if alpha < 0
        eta = alpha/2. + sqrt(alpha*alpha/4.-alpha);
    else
        eta = alpha/2.;
    end
end
    