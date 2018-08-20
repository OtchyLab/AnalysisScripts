function undirected=sortundirected(cell, directed_vector)

counter=1;

for i=1:length(cell.rec)
    directed=0;
    for j=1:length(directed_vector)
            if (i==directed_vector(j) || i>104)
                directed=1;
            end
    end
    if (directed==0)
                undirected.rec(1,counter)=cell.rec(i);
                undirected.motif(counter,:)=cell.motif(i,:);
                undirected.syll(counter,:)=cell.syll(i,:);
                undirected.spikes(counter,:)=cell.spikes(i,:);
                counter=counter+1;
                
            
    end

end