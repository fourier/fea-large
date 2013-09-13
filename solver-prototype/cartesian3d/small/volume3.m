function V = volume3(p1,p2,p3,p4)
% V = volume(p1,p2,p3,p4)
% calculate volume of tethraedra by given 4 nodes in form:
% p1[1..3], p2[1..3], p3[1..3], p4[1..3]
M = [1, p1(1), p1(2), p1(3);
     1, p2(1), p2(2), p2(3);		
     1, p3(1), p3(2), p3(3);		
     1, p4(1), p4(2), p4(3);];

V = det(M)/6.0;