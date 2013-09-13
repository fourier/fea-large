function K = integrate1node(el,F,varargin)
% 1-node linear integration
    w = 1;
    a = 1/4.;
    K = w*F(el,a,a,a,varargin{:})/6.;
end
