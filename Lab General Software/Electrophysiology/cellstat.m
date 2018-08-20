function cellstat(cell, controlIn,controlOut,drugIn,drugOut)

prompt = {'Cell number:','Enter start time (24h):','Enter end time (24h):','Enter depth in nucleus (microns):','Enter the channel #:','Enter age of bird on day of recording(dph):'};
dlg_title = 'Input for cell statistic';
num_lines = 1;
def = {'1','13','14','200','3','70'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
cellFile.cellNo=str2num(answer{1});
cellFile.start=str2num(answer{2});
cellFile.end=str2num(answer{3});
cellFile.depth=str2num(answer{4});
cellFile.channel=str2num(answer{5});
cellFile.dph=str2num(answer{6});

%this calculates the instanteneous firing rates  and the ISI from the cell record
counter=1;
control=find (cell.rec>=controlIn & cell.rec<=controlOut);
drug=find (cell.rec>=drugIn & cell.rec<=drugOut);

for i=1:length(control)
    cell.spikes(control(i),find((cell.spikes(control(i),:)<(cell.motif(control(i),1)-0.05)) | cell.spikes(control(i),:)>cell.motif(control(i),2)))=0;
    raw_raster(control(i),1:size(cell.spikes,2))=cell.spikes(control(i),1:size(cell.spikes,2));
    if (i>1)
        spikes_prev=spikes;
    end
    spikes=nonzeros(raw_raster(control(i),:))';
    for j=1:length(spikes)-1
           
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
 ifr_bins=(2:4:798);
 ifr_hist = histc(ifr,ifr_bins);
 ifr_hist_control = ifr_hist/sum(ifr_hist);% (normalized)
 isi_bins=(0:0.00025:0.1);
 isi_hist_control = histc(isi,isi_bins);
 counter=1;
 for i=1:length(drug)
    cell.spikes(drug(i),find((cell.spikes(drug(i),:)<(cell.motif(drug(i),1)-0.05)) | cell.spikes(drug(i),:)>cell.motif(drug(i),2)))=0;
    raw_raster(drug(i),1:size(cell.spikes,2))=cell.spikes(drug(i),1:size(cell.spikes,2));
    if (i>1)
        spikes_prev=spikes;
    end
    spikes=nonzeros(raw_raster(drug(i),:))';
    for j=1:length(spikes)-1
           
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

ifr_bins=(2:4:798);
ifr_hist = histc(ifr,ifr_bins);
ifr_hist_drug = ifr_hist/sum(ifr_hist);% (normalized)
isi_bins=(0:0.00025:0.1);
isi_hist_drug = histc(isi,isi_bins);

cellFile.ISIDrug{1,1}=ifr_hist_drug;
cellFile.ISIDrug{1,2}=ifr_bins;
cellFile.ISINoDrug{1,1}=ifr_hist_control;
cellFile.ISINoDrug{1,2}=ifr_bins;

time(1:200)=1./ifr_bins;
isiCellsDrugTimeWeighted(1:200)=ifr_hist_drug.*time;
isiCellsNoDrugTimeWeighted(1:200)=ifr_hist_control.*time;
for k=1:200
    activeDrug(k)=sum(isiCellsDrugTimeWeighted(k:200))/sum(isiCellsDrugTimeWeighted(1:200));
    activeNoDrug(k)=sum(isiCellsNoDrugTimeWeighted(k:200))/sum(isiCellsNoDrugTimeWeighted(1:200));
end

cellFile.activeDrug{1,1}=activeDrug(1,:);
cellFile.activeDrug{1,2}=ifr_bins;
cellFile.activeNoDrug{1,1}=activeNoDrug(1,:);
cellFile.activeNoDrug{1,2}=ifr_bins;
cellFile.frDrug{1,1}=[];
cellFile.frDrug{1,2}=[];
cellFile.spontSpikesDrug=[];
cellFile.spontSpikesNoDrug=[];
showCellFile(cellFile);
[FileName,PathName]=uiputfile('*.mat','Save Cell File As:','Statistic_Cell');
save([PathName FileName],'cellFile');
 
