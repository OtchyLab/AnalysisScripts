%%
%Script to processes frequency bands in the aligned songs of pCAF shifted
%birds
clear all;

%Load dataset for day 1
load('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur696\Pur696 120514\Pur696_120514_dataset (platonic,Ch3).mat', 'data')
%load('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Larger Dataset\Pur696\Pur696 120514\Pur696_120514_dataset (platonic, Ch1).mat')
data1 = data;

%Select target syllable bounds from day 1 data
a = transpose(data.templatesyllBreaks);
a = a(:);
buffer = a(1);
band = [a(5),a(6)];

%Load dataset for day 2
load('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Pur696\Pur696 120517\Pur696_120517_dataset (05-14 platonic,Ch3).mat', 'data')
%load('C:\Users\Tim\Desktop\For Ali Paper Test Alignments\Larger Dataset\Pur696\Pur696 120517\Pur696_120517_dataset (platonic, Ch1).mat')
data2 = data;
clear data;

%%
%Go through each rendition and warp the raw audio with the previously found
%path.
rendLength = size(data1.aligned_audioCube,3); 
renditions1 = size(data1.aligned_audioCube,1);
renditions2 = size(data2.aligned_audioCube,1);

%Set up pitch target parameters
pitchCenter = 522; %mean pitch for target syllable
%pitchStd = 12.5; 
numHarm = 10; %number of harmonics to include

%Set up interpolated frequency bins
[~,freqBins,~,~] = spectrogram((data1.PPaudio{1}/(sqrt(mean(data1.PPaudio{1}.^2)))),220,220-44,4096,44150);

%Set up indices for the target pitch harmonics
indx = [];
for  i = 1:numHarm
    diffVect = abs(freqBins-(i*pitchCenter));
    [~,pntr] = min(diffVect);
    indx = [indx;pntr];
end

offset = round(mean(diff(indx)/2));
indx_off = indx+offset;

%%
%Generate hi-res spectrograms and warp using saved paths
intervals_day1 = [];
for i = 1:renditions1
    path = [data1.p{i},data1.q{i}];
    rendLength = size(data1.audioSpecs_LN{i},2);
    
    [warpedOut] = getWarpedStarts(path,a(:));
    intervals_day1 = [intervals_day1;diff(warpedOut)'];
    
    [S,F,T,P] = spectrogram((data1.PPaudio{i}/(sqrt(mean(data1.PPaudio{i}.^2)))),220,220-44,4096,44150);
    
    [alignedSpecAudio1(i,:,:)] = alignSeries(abs(P(1:1000,:)),round(path));
    
    [alignedNeuroPow1(i,:)] = alignRawTS(data1.PPneuro{i},rendLength,round(path));
    
end
%Log-transform and normalize mean spectrogram
trans = 10*log10(squeeze(mean(alignedSpecAudio1,1)));
meanSpecs1 = trans./mean(trans(:));

%Get similarly aligned/warped neural trace
% av_neuroTS1=squeeze(mean(data1.aligned_neuroTS.^2,1)); 
av_neuroTS1=squeeze(mean(alignedNeuroPow1.^2,1)); 
av_windowed_neuroTS1 = mean(running_windows(av_neuroTS1,220,44),2)';
 
intervals_day2 = [];
for i = 1:renditions2
    path = [data2.p{i},data2.q{i}];
    rendLength = size(data2.audioSpecs_LN{i},2);
    
    [warpedOut] = getWarpedStarts(path,a(:));
    intervals_day2 = [intervals_day2;diff(warpedOut)'];
    
    [S,F,T,P] = spectrogram((data2.PPaudio{i}/(sqrt(mean(data2.PPaudio{i}.^2)))),220,220-44,4096,44150);
    
    [alignedSpecAudio2(i,:,:)] = alignSeries(abs(P(1:1000,:)),round(path));
    
    [alignedNeuroPow2(i,:)] = alignRawTS(data2.PPneuro{i},rendLength,round(path));
end
%Log-transform and normalize mean spectrogram
trans = 10*log10(squeeze(mean(alignedSpecAudio2,1)));
meanSpecs2 = trans./mean(trans(:));

%Get similarly aligned/warped neural trace
%av_neuroTS2=squeeze(mean(data2.aligned_neuroTS.^2,1));
av_neuroTS2=squeeze(mean(alignedNeuroPow2.^2,1)); 
av_windowed_neuroTS2 = mean(running_windows(av_neuroTS2,220,44),2)';

%Unwarp spectrograms and neural traces using inverse paths
invPath1 = fliplr([a, round([buffer cumsum(mean(intervals_day1,1))+buffer])']);
invPath2 = fliplr([a, round([buffer cumsum(mean(intervals_day2,1))+buffer])']);

un_meanSpecs1 = alignSeriesSTW(meanSpecs1,invPath1);
un_meanSpecs2 = alignSeriesSTW(meanSpecs2,invPath2);

un_av_neuroTS1 = alignSeriesSTW(av_windowed_neuroTS1,invPath1);
un_av_neuroTS2 = alignSeriesSTW(av_windowed_neuroTS2,invPath2);

%%
%For the template aligned versions
harmonicPowerON1 =  meanSpecs1(indx,:); 
harmonicPowerOFF1 =  meanSpecs1(indx_off,:); 

harmonicPowerON2 =  meanSpecs2(indx,:); 
harmonicPowerOFF2 =  meanSpecs2(indx_off,:); 

temp1 = [];
temp2 = [];
for i=1:numHarm
    temp1(i,:) = harmonicPowerON1(i,:)./harmonicPowerOFF1(i,:);
    temp2(i,:) = harmonicPowerON2(i,:)./harmonicPowerOFF2(i,:);
end

harmonicPower1 = mean(temp1,1);
harmonicPower2 = mean(temp2,1);

%For the unwarped versions
un_harmonicPowerON1 =  un_meanSpecs1(indx,:); 
un_harmonicPowerOFF1 =  un_meanSpecs1(indx_off,:); 

un_harmonicPowerON2 =  un_meanSpecs2(indx,:); 
un_harmonicPowerOFF2 =  un_meanSpecs2(indx_off,:); 

temp1 = [];
temp2 = [];
for i=1:numHarm
    temp1(i,:) = un_harmonicPowerON1(i,:)./un_harmonicPowerOFF1(i,:);
    temp2(i,:) = un_harmonicPowerON2(i,:)./un_harmonicPowerOFF2(i,:);
end

un_harmonicPower1 = mean(temp1,1);
un_harmonicPower2 = mean(temp2,1);

%%
%Plot the template aligned output for review
figure(107);
subplot(6,1,1)
hold on
imagesc(-meanSpecs1(1:560,:)); axis xy
line([band(1),band(1)],[0,560],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
line([band(2),band(2)],[0,560],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
set(gca,'XTick',[],'YTick',[]);
title('Day 0 Spectrogram')

subplot(6,1,2)
hold on
imagesc(-meanSpecs2(1:558,:)); axis xy
line([band(1),band(1)],[0,560],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
line([band(2),band(2)],[0,560],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
set(gca,'XTick',[],'YTick',[]);
title('Day 4 Spectrogram')

subplot(6,1,3)
hold on
plot(smooth(harmonicPower1,5))
plot(smooth(harmonicPower2,5),'r')
line([band(1),band(1)],[.85,1.15],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
line([band(2),band(2)],[.85,1.15],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight;set(gca,'XTick',[]);
title('Harmonic Power Ratio (@ 520Hz)')

subplot(6,1,4)
hold on
diffPlot = (smooth(harmonicPower1,5)-smooth(harmonicPower2,5)).^2;
line([band(1),band(1)],[0,.03],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
line([band(2),band(2)],[0,.03],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
area(diffPlot)
axis tight;set(gca,'XTick',[]);
title('Square Difference in Harmonic Power Ratio')

subplot(6,1,5)
hold on
plot(mat2gray(av_windowed_neuroTS1))
plot(mat2gray(av_windowed_neuroTS2),'r')
line([band(1),band(1)],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
line([band(2),band(2)],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight;set(gca,'XTick',[]);
title('Aligned Neural Traces Before/After pCAF')

subplot(6,1,6)
hold on
diffTrace = (mat2gray(av_windowed_neuroTS1)-mat2gray(av_windowed_neuroTS2)).^2;
area(diffTrace)
axis tight;ylim([0,1])
xlabel('Template Time (ms)')
title('Square Difference of Neural Traces Before/After pCAF')

%%
%Plot the unwarped output
time1 = 1:size(un_meanSpecs1,2);
time2 = 1:size(un_meanSpecs2,2);

shorter = min(length(un_harmonicPower1),length(un_harmonicPower2));
longer = max(length(un_harmonicPower1),length(un_harmonicPower2));

[warpBand1] = interp1(invPath1(:,2),invPath1(:,1),band)-buffer;
[warpBand2] = interp1(invPath2(:,2),invPath2(:,1),band)-buffer;

tShift1 = time1-warpBand1(1); % Absolute time scale for syll-aligned
tShift2 = time2-warpBand2(1)+2; %

minEdge = max(tShift1(1),tShift2(1));
maxEdge = min(tShift1(end),tShift2(end));

%Find the proper starts and stops in the time series
minEdgeIND1 = find(tShift1==minEdge);
minEdgeIND2 = find(tShift2==minEdge);

maxEdgeIND1 = find(tShift1==maxEdge);
maxEdgeIND2 = find(tShift2==maxEdge);

Indx1 = minEdgeIND1:maxEdgeIND1;
Indx2 = minEdgeIND2:maxEdgeIND2;

figure(269);
% subplot(6,1,1)
% hold on
% imagesc(-un_meanSpecs1(1:560,:)); axis xy
% %line([warpBand1(1),warpBand1(1)],[0,size(un_meanSpecs1,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% %line([warpBand1(2),warpBand1(2)],[0,size(un_meanSpecs1,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% axis tight
% %xlim([0,longer])
% set(gca,'XTick',[],'YTick',[]);
% title('Day 0 Spectrogram')

subplot(6,1,1)
hold on
imagesc(-un_meanSpecs2(1:600,:)); axis xy
axis tight
%xlim([-4,longer])
set(gca,'XTick',[],'YTick',[],'TickDir','out');
set(gca,'Box','off')
title('Song Spectrogram')


subplot(6,1,2)
hold on
plot(tShift1(Indx1),smooth(un_harmonicPower1(Indx1),5))
plot(tShift2(Indx2),smooth(un_harmonicPower2(Indx2),5),'r')
line([0,0],[0.85,1.15],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight;set(gca,'XTick',[]);set(gca,'TickDir','out')
set(gca,'Box','off')
title('Harmonic Power Ratio (@ 520Hz)')

subplot(6,1,3)
hold on
diffPlot = (smooth(un_harmonicPower1(Indx1),5)-smooth(un_harmonicPower2(Indx2),5)).^2;
area(tShift1(Indx1),diffPlot)
line([0,0],[0,0.03],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
xlim([0,longer])
axis tight
set(gca,'XTick',[]);
set(gca,'TickDir','out')
set(gca,'Box','off')
title('Square Difference in Harmonic Power Ratio')
 
subplot(6,1,4)
hold on
plot(tShift1(Indx1),mat2gray(un_av_neuroTS1(Indx1)))
plot(tShift2(Indx2),mat2gray(un_av_neuroTS2(Indx2)),'r')
line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% line([band(2),band(2)],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
set(gca,'XTick',[]);
set(gca,'TickDir','out')
set(gca,'Box','off')
title('Aligned Neural Traces Before/After pCAF')
 
subplot(6,1,5)
hold on
diffTrace = (mat2gray(un_av_neuroTS1(Indx1))-mat2gray(un_av_neuroTS2(Indx2))).^2;
area(tShift1(Indx1),diffTrace)
line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
ylim([0,1])
set(gca,'TickDir','out')
set(gca,'Box','off')
xlabel('Absolute Time (ms)')
title('Square Difference Neural Traces Before/After pCAF')

subplot(6,1,6)
hold on
rho = LocalPearson(mat2gray(un_av_neuroTS1(Indx1)),mat2gray(un_av_neuroTS2(Indx2)));
rhoINT = interp1(1:length(rho),rho,linspace(1,length(rho),length(Indx1)));
plot(tShift1(Indx1),rhoINT)
axis tight
ylim([0,1])
set(gca,'YTick',[0,.5,1])
set(gca,'TickDir','out','Box','off')
xlabel('Template Time (ms)')



% diffTrace = (mat2gray(un_av_neuroTS1(Indx1))-mat2gray(un_av_neuroTS2(Indx2))).^2;
% area(tShift1(Indx1),diffTrace)
% line([0,0],[0,1],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% axis tight
% ylim([0,1])
% set(gca,'TickDir','out')
% set(gca,'Box','off')
% xlabel('Absolute Time (ms)')
% title('Square Difference Neural Traces Before/After pCAF')




%%
%Plot the unwarped, syllable sligned data
longer = max(size(un_meanSpecs1(:,warpBand1(1):end),2),size(un_meanSpecs2(:,warpBand2(1):end),2));

figure(300);
subplot(4,1,1)
hold on
imagesc(-un_meanSpecs1(:,warpBand1(1):end)); axis xy
%[warpBand1] = interp1(invPath1(:,2),invPath1(:,1),band)-buffer;
%line([warpBand1(1),warpBand1(1)],[0,size(un_meanSpecs1,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
%line([warpBand1(2),warpBand1(2)],[0,size(un_meanSpecs1,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
xlim([0,longer])
set(gca,'XTick',[],'YTick',[]);
title('Day 0 Spectrogram')

subplot(4,1,2)
hold on
imagesc(-un_meanSpecs2(:,warpBand2(1):end)); axis xy
%[warpBand2] = interp1(invPath2(:,2),invPath2(:,1),band)-buffer;
%line([warpBand2(1),warpBand2(1)],[0,size(un_meanSpecs2,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
%line([warpBand2(2),warpBand2(2)],[0,size(un_meanSpecs2,1)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight
xlim([0,longer])
set(gca,'XTick',[],'YTick',[]);
title('Day 4 Spectrogram')

subplot(4,1,3)
hold on
plot(smooth(un_harmonicPower1(warpBand1(1):end),5))
plot(smooth(un_harmonicPower2(warpBand2(1):end),5),'r')
% line([band(1),band(1)],[.85,1.15],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% line([band(2),band(2)],[.85,1.15],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
axis tight;set(gca,'XTick',[]);
ylim([0.85,1.15])
title('Harmonic Power Ratio (@ 520Hz)')

subplot(4,1,4)
hold on
shorter = min(length(un_harmonicPower1(warpBand1(1):end)),length(un_harmonicPower2(warpBand2(1):end)));
longer = max(length(un_harmonicPower1(warpBand1(1):end)),length(un_harmonicPower2(warpBand2(1):end)));
diffPlot = (smooth(un_harmonicPower1(warpBand1(1):warpBand1(1)+shorter-1),5)-smooth(un_harmonicPower2(warpBand2(1):warpBand2(1)+shorter-1),5)).^2;
% line([band(1),band(1)],[0,.03],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
% line([band(2),band(2)],[0,.03],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
plot(diffPlot)
axis tight
xlim([0,longer])
ylim([0,0.03])
xlabel('Absolute Time (ms)')
title('Power Difference Squared')