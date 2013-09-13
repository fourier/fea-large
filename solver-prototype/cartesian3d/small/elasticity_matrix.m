function M = elasticity_matrix(F)
    M = elasticity_A5();
end

function M = elasticity_linear()
% function M = elasticity_matrix()
% Returns Elasticity matrix 6x6 for
% the model A5
% A5: lambda*I1(Epsilon)+2mu*Epsilon

% turn on for elastic model
%if not(true)
    E = 1e9;
    nu = 0.3;
    mul = E/((1+nu)*(1-2*nu));
    M = mul*[1-nu, nu, nu, 0, 0, 0;
             nu, 1-nu, nu, 0, 0, 0;
             nu, nu, 1-nu, 0, 0, 0;
             0,	0,	0,  (1-2*nu)/2., 0, 0;
             0, 0,  0,  0, (1-2*nu)/2., 0;
             0, 0,  0,  0, 0, (1-2*nu)/2.];

end


function M = elasticity_A5()
% model A5
    l = 100;
    m = 100;
    M = [l+2*m, l, l, 0, 0, 0;
         l,	l+2*m, l, 0, 0, 0;
         l, l, l+2*m, 0, 0, 0;
         0,	0,	0,	2*m, 0, 0;
         0, 0,  0,  0, 2*m, 0;
         0, 0,  0,  0, 0, 2*m];
end

function M = elasticity_A5_large(F)
    M = zeros(6,6);
    
end