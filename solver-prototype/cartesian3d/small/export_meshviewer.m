function export_meshviewer(filename,nodes,elements,IDX,X) 

f = fopen(filename,'w+');
str = sprintf('%d 3 0\n', size(nodes,1));
fprintf(f,str);
for i = 1:size(nodes,1)
  str = sprintf('%d %f %f %f\n', i,nodes(i,1),nodes(i,2),nodes(i,3));
  fprintf(f,str);
end
str = sprintf('%d 4 6 Sxx Syy Szz Sxy Sxz Syz\n', size(elements,1));
fprintf(f,str);
for i = 1:size(elements,1)
  s = stress(nodes,elements,IDX,i,X);
  str = sprintf('%d %d %d %d %d %f %f %f %f %f %f\n',...
    i,elements(i,1),elements(i,2),elements(i,3),elements(i,4),...
    s(1),s(2),s(3),s(4),s(5),s(6));
  fprintf(f,str);
end
fclose(f);