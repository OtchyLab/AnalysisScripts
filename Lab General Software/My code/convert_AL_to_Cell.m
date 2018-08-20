function totalCell=convert_AL_to_Cell(idata)
counter(1:200)=1;   %motif counter 
for i=1:137
    if (length(idata(1,i).cells)>0)
        for j=1:idata(1,i).sylnum-3
            if (strcmp(idata(1,i).sname(j:j+3),'abcd'))
                for cellnum=1:length(idata(1,i).cells)
                    spikes=[];
                    totalCell.rec(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)))=counter(idata(1,i).cells(cellnum));
                    totalCell.motif(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),1)=idata(1,i).sindex(j,1);
                    totalCell.motif(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),2)=idata(1,i).sindex(j+3,2);
                    for syllno=1:4
                        totalCell.syll(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),2*syllno-1)=idata(1,i).sindex(j+syllno-1,1);
                        totalCell.syll(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),2*syllno)=idata(1,i).sindex(j+syllno-1,2);
                    end
                    if (cellnum==1)
                        spikes(1:idata(1,i).nspikes(cellnum))=idata(1,i).spikes(1:idata(1,i).nspikes(1));
                    else
                        spikes(1:idata(1,i).nspikes(cellnum)-idata(1,i).nspikes(cellnum-1)+1)=idata(1,i).spikes(idata(1,i).nspikes(cellnum-1):idata(1,i).nspikes(cellnum));
                    end
                    if (totalCell.motif(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),1)-0.3*40000<0)
                        startmotif=1;
                    else
                        startmotif=find(spikes>floor(totalCell.motif(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),1)-0.3*40000),1,'first');
                    end
                    endmotif=find(spikes<floor(totalCell.motif(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),2)+0.3*40000),1,'last');
                    totalCell.spikes(idata(1,i).cells(cellnum),counter(idata(1,i).cells(cellnum)),1:endmotif-startmotif+1)=spikes(startmotif:endmotif);
                    counter(idata(1,i).cells(cellnum))=counter(idata(1,i).cells(cellnum))+1;
                end
            end
        end
    else
        %do nothing
    end
end
   for i=1:55
            presentCell=[];
            endfile=find(totalCell.rec(i,:)>0,1,'last');
            presentCell.rec(1:endfile)=totalCell.rec(i,1:endfile);
            presentCell.motif(1:endfile,1:2)=totalCell.motif(i,1:endfile,1:2)/40000;
            presentCell.syll(1:endfile,1:8)=totalCell.syll(i,1:endfile,1:8)/40000;
            presentCell.spikes(1:endfile,:)=totalCell.spikes(i,1:endfile,:)/40000;
            
            save(['ALCell' num2str(i)],'presentCell');
   end
                    
                    
                    