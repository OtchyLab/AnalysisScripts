function [ifr_hist,ifr_bins,isi_hist,isi_bins]=ISI_IFR_cell(cell, filename)

%this calculates the instanteneous firing rates  and the ISI from the cell record
counter=1;
for i=1:size(cell.spikes,1)
    cell.spikes(i,find((cell.spikes(i,:)<(cell.motif(i,1)-0.05)) | cell.spikes(i,:)>cell.motif(i,2)))=0;
    raw_raster(i,1:size(cell.spikes,2))=cell.spikes(i,1:size(cell.spikes,2));
    if (i>1)
        spikes_prev=spikes;
    end
    spikes=nonzeros(raw_raster(i,:))';
    for (j=1:length(spikes)-1)
           
        if (i>1) 
            testx=find(spikes_prev == spikes(j)); %test if the spike in the motif was already 
            %accounted for in the previous motif (this since the cell array takes spikes 300 ms 
            %before and after a motif)
            if (testx>0)
                %do nothing: this is a duplicate spike from previous motif
            else
                spikes(j)=spikes(j)+rand*(1/40000); %rand is there to make the distribution smoother
                ifr(counter)=1/(spikes(j+1)-spikes(j));
                isi(counter)=spikes(j+1)-spikes(j);
                counter=counter+1;
            
            end
        end
    end
end
 ifr_bins=(0:4:800);
 ifr_hist = histc(ifr,ifr_bins);
 ifr_hist = ifr_hist/sum(ifr_hist);% (normalized)
 figure;plot(ifr_bins,ifr_hist);
 isi_bins=(0:0.00025:0.1);
 isi_hist = histc(isi,isi_bins);
 figure;semilogy(isi_bins,isi_hist);
 save (filename,'isi_bins','isi_hist','ifr_bins','ifr_hist');
end
