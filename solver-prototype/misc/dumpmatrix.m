function dumpmatrix(filename,m)
f = fopen(filename,'w+');
for i = 1:size(m,1)
	for j = 1:size(m,2)
		str = sprintf('%.8f ',m(i,j));
		fprintf(f,str);
	end
	fprintf(f,'\n');
end
fclose(f);
