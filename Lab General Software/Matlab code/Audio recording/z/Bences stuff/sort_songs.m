%sort.m  This sorts the songs in the present directory. Assumes that
%"Analyze" is running
    i=1;
    j=1;
    clear junk_vector song_vector song_hist junk_hist
for k=1:size(namehelpsA,2) %go over all the files in the present directory

    filenameA=[loadpath birdnameA '/' dirnameA '/' namehelpsA{k}];
    yAInfo=daqread(filenameA,'info');
    yA=daqread(filenameA)';
    
    [spectra,f,t]=specgram1(yA(1,:),512,scanrateA,400,200); ylim([0 12000]);
    spectra=abs(spectra);
    ratio=PSRatio(scanrateA,spectra,t,f);
    
    %is this a song??
    time=yAInfo.ObjInfo.InitialTriggerTime(4)*60+yAInfo.ObjInfo.InitialTriggerTime(5); %time in minutes
    tcrossings;
    if (song)
       song_vector(i)=time;
       i=i+1;
    else
       junk_vector(j)=time;
       j=j+1;
    end
 
end
time_int= 0:1:1440;
song_hist = hist(song_vector,time_int);
junk_hist = hist(junk_vector,time_int);
figure; subplot(2,1,1);plot(song_hist,'color','red'); subplot(2,1,2); plot(junk_hist,'color','blue'); 