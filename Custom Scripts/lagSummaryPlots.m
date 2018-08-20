function lagSummaryPlots

%Get the Processed file to analyze
[temp,path] = uigetfile('*.mat','Select the processed data file.','MultiSelect','on');

if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames.name = temp;
end

%Define variables

maxLag = -1;
snipLags = -100:100;
%Loop through each day's data
for i = 1:length(fnames)
    %Load the data
    load([path filesep fnames(i).name],'m_corrSnip','corrFull','lagsFull')    
    %fullCorr{i} = mat2gray(corrFull);
    %fullCorr{i} = (corrFull);
    fullCorr{i} = corrFull/max(corrFull);
    fullLags{i} = lagsFull-7;
    %snipCorr(i,:) = mat2gray(m_corrSnip);
    %snipCorr(i,:) = (m_corrSnip);
    snipCorr(i,:) = m_corrSnip/max(m_corrSnip);
    
    maxLag = max(maxLag,length(fullLags{i}));
end
m_snipCorr = mean(snipCorr(i,:),1);

%Align full corr arrays for averaging
corrBlock = NaN(i,maxLag);
lagsBlock = NaN(i,maxLag);
for i = 1:length(fnames)
    pntr = ((maxLag-length(fullLags{i}))/2);
    lagsBlock(i,pntr+1:(maxLag-pntr)) = fullLags{i};
    %corrBlock(i,pntr+1:(end-pntr)) = fullCorr{i};
    corrBlock(i,pntr+1:(end-pntr)) = fullCorr{i};
end

m_lags = nanmean(lagsBlock,1);
m_corr = nanmean(corrBlock,1);

%Plot the outputs
fig1 = figure;
subplot(2,1,1)
plot(m_lags,m_corr)
title('Correlation Lags Between Song and Neural Activity')
xlabel('Neural Signal Lead (ms)')
ylabel('Correlation')
set(gca,'TickDir','out','Box','off')
%ylim([0,1])

subplot(2,1,2)
plot(snipLags,m_snipCorr)
xlabel('Neural Signal Lead (ms)')
ylabel('Correlation')
set(gca,'TickDir','out','Box','off')
%ylim([0,1])

%Save outputs to file
save([path fnames(i).name(1:6) ' over days peak normalized.mat'],'snipLags','m_snipCorr','m_lags','m_corr')
saveas(fig1,[path fnames(i).name(1:6) ' over days peak normalized plot.fig'])



