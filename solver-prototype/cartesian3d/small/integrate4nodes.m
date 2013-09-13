function K = integrate4nodes(el,F,varargin)
% 4-nodes quadratic integration
    a = 0.58541020;
    b = 0.13819660;
    w = 1/4.;
    K = w*(F(el,a,b,b,varargin{:})+F(el,b,a,b,varargin{:})+ ...
           F(el,b,b,a,varargin{:})+F(el,b,b,b,varargin{:}))/6.;
end
