function h = test_isoform(r,s,t)
h = 0;
for i = 1:10
	h = h + isoform(i,r,s,t);
end