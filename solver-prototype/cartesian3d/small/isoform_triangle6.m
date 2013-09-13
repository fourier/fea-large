function h = isoform_triangle6(i,r,s)
% function h = isoform_triangle6(i,r,s)
% function for calculation value of shape function for 6-noded
% triangle
% by given node number i and local coordinates r,s 
% where r,s from [0;1]
% all functions are taken from the book:
% "The Finite Element Method for 3D Thermomechanical Applications - Guido Dhond" p.72

  h = 0;
  if i == 1
    h = (2*(1-r-s-t)-1)*(1-r-s-t);
  elseif i == 2
    h = (2*r-1)*r;
  elseif i == 3
    h = (2*s-1)*s;
  elseif i == 4
    h = (2*t-1)*t;
  elseif i == 5
    h = 4*r*(1-r-s-t);
  elseif i == 6
    h = 4*r*s;
  elseif i == 7
    h = 4*s*(1-r-s-t);
  elseif i == 8
    h = 4*t*(1-r-s-t);
  elseif i == 9
    h = 4*r*t;
  else % i == 10
    h = 4*s*t;
  end

