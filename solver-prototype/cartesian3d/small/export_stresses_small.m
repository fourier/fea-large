function export_stresses_small(filename,nodes0,elements,IDX,U)
% function export_stresses_small(filename,nodes0,elements,IDX,U)
    file = fopen(filename,'w+');
    fprintf(file,'N\tSxx\tSyy\tSzz\tSxy\tSxz\tSyz\n');
    for i = 1:size(elements,1)
        q = create_q_from_u(IDX,U,i);
        el = element(nodes0,elements,i);
        a = 0.58541020;
        b = 0.13819660;
        r = a;
        s = b;
        t = b;
        e = epsilon_lin(el,q,r,s,t);
        A = elasticity_matrix();
        s = A*e;
        fprintf(file,'%d\t%e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',...
                i,s(1),s(2),s(3),s(4),s(5),s(6));
    end 
    fclose(file);
end