function MUA_Stopping
clear all

%Request and load dataset
[fname,pathLoc] = uigetfile('C:\Users\Tim\Desktop\Stopping Analysis\*.mat');
tots = [pathLoc, fname];
load(tots,'data')

%Create indices for each of the analysis groups
compMotif = find(data.postSyll==5);
stopMotif = find(data.postSyll==99);

%Get the important syllable breaks
a = transpose(data.templatesyllBreaks);
a = a(:);

%Interval Durations from the paths and template files
rendNum = length(data.audio);
IntMat = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path = [data.p{i},data.q{i}];
    
    %Get syllable breaks for the rendition
    sylBreaks(i,:) = getWarpedStarts(path,a(:));
    %IntMat = [IntMat; diff(warpedOut)'];
    
    %Gather unwarped audio and neuro power
    un_audio{i} = -.1*sum(data.audioSpecs_LN{i},1);
    %un_neuro{i} = mean(running_windows(data.PPneuro{i}.^2,220,44),2);
    smd = smooth(data.PPneuro{i}.^2,220,'moving')';
    un_neuro{i} = interp1(1:length(smd),smd,linspace(1,length(smd),length(un_audio{i})));
end

%Plot output
fig1 = figure;
tarLength = (a(6)-a(5));
maxLen = 0;
for i=compMotif
    rendLength = sylBreaks(i,6)-sylBreaks(i,5);
    scaleX(i) = tarLength/rendLength;
    
    audio_snip{i} = un_audio{i}(sylBreaks(i,5):end);
    neuro_snip{i} = un_neuro{i}(sylBreaks(i,5):end);
    snipLength = length(audio_snip{i});
    
    audio_Rsnip{i} = interp1(1:snipLength,audio_snip{i},linspace(1,snipLength,round(snipLength*scaleX(i))));
    neuro_Rsnip{i} = interp1(1:snipLength,neuro_snip{i},linspace(1,snipLength,round(snipLength*scaleX(i))));
    xs = (1:length(audio_Rsnip{i}))-tarLength;
    maxLen = max(length(xs),maxLen);
    
    %Audio power traces
    subplot(2,1,1)
    hold on
    plot(xs,audio_Rsnip{i},':b','LineWidth',0.25);
    
    %Neural power traces
    subplot(2,1,2)
    hold on
    plot(xs,neuro_Rsnip{i},':b','LineWidth',0.25);
end

for i=stopMotif
    rendLength = sylBreaks(i,6)-sylBreaks(i,5);
    scaleX(i) = tarLength/rendLength;
    
    audio_snip{i} = un_audio{i}(sylBreaks(i,5):end);
    neuro_snip{i} = un_neuro{i}(sylBreaks(i,5):end);
    snipLength = length(audio_snip{i});
    
    audio_Rsnip{i} = interp1(1:snipLength,audio_snip{i},linspace(1,snipLength,round(snipLength*scaleX(i))));
    neuro_Rsnip{i} = interp1(1:snipLength,neuro_snip{i},linspace(1,snipLength,round(snipLength*scaleX(i))));
    xs = (1:length(audio_Rsnip{i}))-tarLength;
    maxLen = max(length(xs),maxLen);
    
    %Audio power traces
    subplot(2,1,1)
    hold on
    plot(xs,audio_Rsnip{i},':r','LineWidth',0.25);
    
    %Neural power traces
    subplot(2,1,2)
    hold on
    plot(xs,neuro_Rsnip{i},':r','LineWidth',0.25);
end

audioMat = NaN(rendNum,maxLen);
neuroMat = NaN(rendNum,maxLen);
for i=1:rendNum
    %Reformat for stats
    audioMat(i,1:length(audio_Rsnip{i})) = audio_Rsnip{i};
    neuroMat(i,1:length(neuro_Rsnip{i})) = neuro_Rsnip{i};
end
xs = (1:maxLen)-tarLength;

%Audio power traces
subplot(2,1,1)
hold on
plot(xs,nanmean(audioMat(compMotif,:),1),'k','LineWidth',2);
plot(xs,nanmean(audioMat(stopMotif,:),1),'g','LineWidth',2);

%Neural power traces
subplot(2,1,2)
hold on
plot(xs,nanmean(neuroMat(compMotif,:),1),'k','LineWidth',2);
plot(xs,nanmean(neuroMat(stopMotif,:),1),'g','LineWidth',2);

%Format Plots
subplot(2,1,1)
line([0,0],[-20,-80],'Color','k','LineStyle','--')
set(gca,'XTick',[],'XTickLabel',[])
set(gca,'TickDir','out','Box','off')
ylabel('Audio Power (dB)')
tits = ['Unwarped Aligned Audio and Neuro Power for ' fname(1:end-11) ' ' fname((end-7):(end-4))];
title(tits,'Interpreter','none')

subplot(2,1,2)
line([0,0],[0,0.5],'Color','k','LineStyle','--')
set(gca,'TickDir','out','Box','off')
ylim([0,0.5])
xlabel('Time (ms)')
ylabel('Neuro Power (V^2)')

%Save data to file
clear('data')
savename = [pathLoc 'Aligned\' fname(1:end-11), fname((end-7):(end-4)) ' alignedVar.mat'];
save(savename,'*')

saveas(fig1,[pathLoc 'Aligned\' fname(1:end-11), fname((end-7):(end-4)) ' aligned.fig'])
saveas(fig1,[pathLoc 'Aligned\' fname(1:end-11), fname((end-7):(end-4)) ' aligned.eps'])


