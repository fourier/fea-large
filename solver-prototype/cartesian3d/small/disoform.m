function h = disoform(shape,dof,r,s,t)
% function h = disoform(i,r,s,t)
% function for calculation derivatives of shape 
% function of 10noded tetrahedra element
% with respect to local coordinate system
% shape - number of node(and corresponding shape function)
% dof - degree of freedom, dof = 1 is r, dof = 2 is s, dof = 3 is t
% r,s,t is [0;1] - local coordinates

  h = 0;
  if dof == 1 % r
    h = df_dr(shape,r,s,t);
  elseif dof == 2 % s
    h = df_ds(shape,r,s,t);
  else % dof == 3, t
    h = df_dt(shape,r,s,t);
  end

% 
% Particular derivatives
%

function h = df_dr(i,r,s,t)
  h = 0;
  if i == 1
    h = 4*t+4*s+4*r-3;
  elseif i == 2
    h = 4*r-1;
  elseif i == 3
    h = 0;
  elseif i == 4
    h = 0;
  elseif i == 5
    h = -4*t-4*s-8*r+4;
  elseif i == 6
    h = 4*s;
  elseif i == 7
    h = -4*s;
  elseif i == 8
    h = -4*t;
  elseif i == 9
    h = 4*t;
  else % i == 10
    h = 0;
  end

function h = df_ds(i,r,s,t)
  h = 0;
  if i == 1
    h = 4*t+4*s+4*r-3;
  elseif i == 2
    h = 0;
  elseif i == 3
    h = 4*s-1;
  elseif i == 4
    h = 0;
  elseif i == 5
    h = -4*r;
  elseif i == 6
    h = 4*r;
  elseif i == 7
    h = -4*t-8*s-4*r+4;
  elseif i == 8
    h = -4*t;
  elseif i == 9
    h = 0;
  else % i == 10
    h = 4*t;
  end

function h = df_dt(i,r,s,t)
  h = 0;
  if i == 1
    h = 4*t+4*s+4*r-3;
  elseif i == 2
    h = 0;
  elseif i == 3
    h = 0;
  elseif i == 4
    h = 4*t-1;
  elseif i == 5
    h = -4*r;
  elseif i == 6
    h = 0;
  elseif i == 7
    h = -4*s;
  elseif i == 8
    h = -8*t-4*s-4*r+4;
  elseif i == 9
    h = 4*r;
  else % i == 10
    h = 4*s;
  end

