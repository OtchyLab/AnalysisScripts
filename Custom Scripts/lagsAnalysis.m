function lagsAnalysis
%This function calulates the lag between audio and neural signals using two
%different methods. The first looks at the whole motif; the second focuses
%just on the 50ms before and after each syllable onset and offset.
%
%The user is first instructed to select the dataset to work from. It is
%important that:
%   (i) This dataset be created (by StretchEm) with 0ms lag
%  (ii) The motif buffer be more than 50ms
%
%The function prints the outputs to screen and saves to file.

%Get the Processed file to analyze
[fname,path] = uigetfile('*.mat','Select the processed data file.');
load([path fname], 'data')

%Calculate the average time series
mAudio = mean(data.audioPower,1);
mNeuro = mean(data.neuroPower,1);

%Get Syllable onset times
sylStarts = data.templatesyllBreaks(:,1)';

%Calculate correlation for the full timeseries
%[corrFull,lagsFull]=xcorr(mat2gray(mAudio),mat2gray(mNeuro),'coeff');
[corrFull,lagsFull]=xcorr(mat2gray(mAudio)-mean(mat2gray(mAudio)),mat2gray(mNeuro)-mean(mat2gray(mNeuro)),'none');
[~,indx] =  max(corrFull);
maxFull = lagsFull(indx);

%Calculate correlation for 100ms Syll-Snips

for i = 1:length(sylStarts)
   %Get the required snips
   snipRange = (sylStarts(i)-50):(sylStarts(i)+50);
   snipAudio(i,:) = mAudio(snipRange);
   snipNeuro(i,:) = mNeuro(snipRange);
   
   %Calculate the snip correlation
   %[corrSnip(i,:),lagsSnip(i,:)]=xcorr(mat2gray(snipAudio(i,:)),mat2gray(snipNeuro(i,:)),'coeff');
   [corrSnip(i,:),lagsSnip(i,:)]=xcorr(mat2gray(snipAudio(i,:))-mean(mat2gray(snipAudio(i,:))),mat2gray(snipNeuro(i,:))-mean(mat2gray(snipNeuro(i,:))),'none');
end
m_corrSnip = mean(corrSnip,1);
[~,indx] =  max(m_corrSnip);
maxSnip = lagsSnip(i,(indx));

%Plot the outputs
fig1 = figure;
subplot(2,1,1)
plot(lagsFull,corrFull)
title('Correlation Lags Between Song and Neural Activity')
xlabel('Neural Signal Lead (ms)')
ylabel('Correlation')
set(gca,'TickDir','out','Box','off')
%ylim([0,1])

subplot(2,1,2)
hold on
for i = 1:length(sylStarts)
    plot(lagsSnip(i,:),corrSnip(i,:),':')
end
plot(lagsSnip(i,:),m_corrSnip,'r')
xlabel('Neural Signal Lead (ms)')
ylabel('Correlation')
set(gca,'TickDir','out','Box','off')
%ylim([0,1])

%Save important stuff
% save([path fname(1:end-4) ' lag summary.mat'],'mAudio','mNeuro','lagsFull','corrFull','maxFull','snipAudio','snipNeuro','corrSnip','lagsSnip','m_corrSnip','maxSnip')
% saveas(fig1,[path fname(1:end-4) ' lag plot.fig'])



