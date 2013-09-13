function dumpvector(filename,v)
f = fopen(filename,'w+');
for i = 1:size(v,1)
	str = sprintf('%.8f\n',v(i));
	fprintf(f,str);
end
fclose(f);