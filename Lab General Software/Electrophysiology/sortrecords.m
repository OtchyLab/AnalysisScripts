function newcell=sortrecords(cell)

[r,index]=sort(cell.rec);
newcell.rec=r;
for i=1:length(cell.rec)
    
    newcell.syll(i,:)=cell.syll(index(i),:);
    newcell.motif(i,:)=cell.motif(index(i),:);
    newcell.spikes(i,:)=cell.spikes(index(i),:);
    
end