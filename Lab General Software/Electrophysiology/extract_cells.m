function cell=extract_cells(clusters,cell_no);

j=1;
for (i=1:length(clusters))
    if (clusters(i,1)==cell_no)
        cell(j)=clusters(i,2);
        j=j+1;
    end
end