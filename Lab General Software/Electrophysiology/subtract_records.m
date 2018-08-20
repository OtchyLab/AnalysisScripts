function newcellfile=subtract_records(cellfile,from,to)
%from and to are inclusive records. these will also be deleted
newcellfile.spikes(1:from-1,:)=cellfile.spikes(1:(from-1),:);
newcellfile.spikes(from:size(cellfile.spikes,1)-(to-from+1),:)=cellfile.spikes((to+1):end,:);
newcellfile.syll(1:from-1,:)=cellfile.syll(1:(from-1),:);
newcellfile.syll(from:size(cellfile.syll,1)-(to-from+1),:)=cellfile.syll((to+1):end,:);
newcellfile.motif(1:from-1,:)=cellfile.motif(1:(from-1),:);
newcellfile.motif(from:size(cellfile.motif,1)-(to-from+1),:)=cellfile.motif((to+1):end,:);
newcellfile.rec=[cellfile.rec(1:from-1) cellfile.rec(to+1:end)];
end