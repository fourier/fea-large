function export_gmsh(undef_filename,def_filename,nodes,elements,IDX,X)
% function export_gmsh(filename,nodes,elements,IDX,X)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Export undeformed file with displacements
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  undef = fopen(undef_filename,'w+');
  write_header(undef);
  write_geometry(undef,nodes,elements);
  write_displacements(undef,nodes,X);
  fclose(undef);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Export undeformed file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  def_nodes = update_with_solution(nodes,X);
  def = fopen(def_filename,'w+');
  write_header(def);
  write_geometry(def,def_nodes,elements);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Postprocessing section
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 
  %% Elements Data postprocessing - section
  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Prepare stresses array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  stresses = zeros(size(elements,1),6);
%  for i = 1:size(elements,1)
%    s = stress(nodes,elements,IDX,i,X);
%    stresses(i,:) = s';
%  end  
%  write_elements_section(f,'Sxx',stresses(:,1));
%  write_elements_section(f,'Syy',stresses(:,2));
%  write_elements_section(f,'Szz',stresses(:,3));
%  write_elements_section(f,'Sxy',stresses(:,4));
%  write_elements_section(f,'Sxz',stresses(:,5));
%  write_elements_section(f,'Syz',stresses(:,6));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 
  %% Element-Nodes Data postprocessing - section
  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Prepare strains/stresses arrays
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nelems = size(elements,1);
  nnodes = size(elements,2);
  strains = zeros(6,nelems*nnodes);
  stresses = zeros(6,nelems*nnodes);
  for i = 1:nelems
    [e,s] = stress_in_nodes(nodes,elements,IDX,i,X);
    strains(:,i*nnodes-nnodes+1:i*nnodes) = e;
    stresses(:,i*nnodes-nnodes+1:i*nnodes) = s;
  end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Export strains array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  write_elementnodes_section(def,'Strain_xx',nnodes,strains(1,:));
  write_elementnodes_section(def,'Strain_yy',nnodes,strains(2,:));
  write_elementnodes_section(def,'Strain_zz',nnodes,strains(3,:));
  write_elementnodes_section(def,'Strain_xy',nnodes,strains(4,:));
  write_elementnodes_section(def,'Strain_xz',nnodes,strains(5,:));
  write_elementnodes_section(def,'Strain_yz',nnodes,strains(6,:));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Export stresses array
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  write_elementnodes_section(def,'Stress_xx',nnodes,stresses(1,:));
  write_elementnodes_section(def,'Stress_yy',nnodes,stresses(2,:));
  write_elementnodes_section(def,'Stress_zz',nnodes,stresses(3,:));
  write_elementnodes_section(def,'Stress_xy',nnodes,stresses(4,:));
  write_elementnodes_section(def,'Stress_xz',nnodes,stresses(5,:));
  write_elementnodes_section(def,'Stress_yz',nnodes,stresses(6,:));
  %%
  %% Close file
  %%
  fclose(def);
end

%%
%% MeshFormat section
%%
function write_header(f)
  fprintf(f,'$MeshFormat\n');
  fprintf(f,'2.0 0 8\n');
  fprintf(f,'$EndMeshFormat\n');
end


function write_geometry(f,nodes,elements)
  %%
  %% Nodes section
  %%
  fprintf(f,'$Nodes\n');
  str = sprintf('%d\n',size(nodes,1));
  fprintf(f,str);
  for i = 1:size(nodes,1)
    str = sprintf('%d %f %f %f\n', i,nodes(i,1),nodes(i,2),nodes(i,3));
    fprintf(f,str);
  end
  fprintf(f,'$EndNodes\n');
  %% 
  %% Elements section
  %%
  fprintf(f,'$Elements\n');
  str = sprintf('%d\n', size(elements,1));
  fprintf(f,str);
  for i = 1:size(elements,1)
    str = sprintf('%d 11 3 1 1 1 %d %d %d %d %d %d %d %d %d %d\n',...
      i,elements(i,1),elements(i,2),elements(i,3),elements(i,4),...
      elements(i,5),elements(i,6),elements(i,7),elements(i,8),...
      elements(i,10),elements(i,9));
    fprintf(f,str);
  end
  fprintf(f,'$EndElements\n');
end

function write_displacements(f,nodes,X)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Nodes Data postprocessing - section
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(f,'$NodeData\n');
  %
  % number-of-string-tags
  % < "string-tag" >
  % ...
  % number-of-real-tags
  % < real-tag >
  % ...
  % number-of-integer-tags
  % < integer-tag >
  % ...
  % elm-number value ...
  % ...
  fprintf(f,'1\n');
  str = '"Ux,Uy,Uz"';
  fprintf(f,str);
  fprintf(f,'1\n');
  fprintf(f,'0\n');
  fprintf(f,'3\n');
  fprintf(f,'0\n');
  fprintf(f,'3\n');
  str = sprintf('%d\n',length(nodes));
  fprintf(f,str);
  for i = 1:length(nodes)
    str = sprintf('%d %f %f %f\n',i,X(i*3-2),X(i*3-1),X(i*3)); 
    fprintf(f,str);
  end  
  fprintf(f,'$EndNodeData\n');
end

function write_elements_section(f,name,values)
  fprintf(f,'$ElementData\n');
  %
  % number-of-string-tags
  % < "string-tag" >
  % ...
  % number-of-real-tags
  % < real-tag >
  % ...
  % number-of-integer-tags
  % < integer-tag >
  % ...
  % elm-number value ...
  % ...
  fprintf(f,'1\n');
  str = sprintf('"%s"\n',name);
  fprintf(f,str);
  fprintf(f,'1\n');
  fprintf(f,'0\n');
  fprintf(f,'3\n');
  fprintf(f,'0\n');
  fprintf(f,'1\n');
  str = sprintf('%d\n',length(values));
  fprintf(f,str);
  for i = 1:length(values)
    str = sprintf('%d %f\n',i,values(i)); 
    fprintf(f,str);
  end  
  fprintf(f,'$EndElementData\n');
end

function write_elementnodes_section(f,name,nnodes,values)
  fprintf(f,'$ElementNodeData\n');
  %
  % number-of-string-tags
  % < "string-tag" >
  % ...
  % number-of-real-tags
  % < real-tag >
  % ...
  % number-of-integer-tags
  % < integer-tag >
  % ...
  % elm-number value ...
  % ...
  fprintf(f,'1\n');
  str = sprintf('"%s"\n',name);
  fprintf(f,str);
  fprintf(f,'1\n');
  fprintf(f,'0\n');
  fprintf(f,'3\n');
  fprintf(f,'0\n');
  fprintf(f,'1\n');
  nelems = length(values)/nnodes;
  str = sprintf('%d\n',nelems);
  fprintf(f,str);
  for i = 1:nelems
    pos = (i-1)*nnodes+1;
    str = sprintf('%d %d %f %f %f %f %f %f %f %f %f %f\n',i,nnodes,...
    values(pos),values(pos+1),values(pos+2),values(pos+3),values(pos+4),...
    values(pos+5),values(pos+6),values(pos+7),values(pos+8),values(pos+9));
    fprintf(f,str);
  end  
  fprintf(f,'$EndElementNodeData\n');
end
