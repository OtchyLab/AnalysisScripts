function lagSummaryOverBirdsPlots

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
    load([path filesep fnames(i).name],'m_corr','m_lags','m_snipCorr')    
    %fullCorr{i} = mat2gray(m_corr);
    fullCorr{i} = m_corr;
    fullLags{i} = m_lags;
    %snipCorr(i,:) = mat2gray(m_snipCorr);
    snipCorr(i,:) = m_snipCorr;
    
    maxSnip = find(m_snipCorr==max(m_snipCorr));
    peakSnip(i) = snipLags(maxSnip);
    
    maxCorr = find(m_corr==max(m_corr));
    peakTime(i) = m_lags(maxCorr);
    
    maxLag = max(maxLag,length(fullLags{i}));
end
m_snipCorr = mean(snipCorr(i,:),1);
sem_snipCorr = std(snipCorr(i,:),1)./sqrt(i);

%Align full corr arrays for averaging
corrBlock = NaN(i,maxLag);
lagsBlock = NaN(i,maxLag);
for i = 1:length(fnames)
    pntr = ((maxLag-length(fullLags{i}))/2);
    lagsBlock(i,pntr+1:(maxLag-pntr)) = fullLags{i};
    corrBlock(i,pntr+1:(end-pntr)) = fullCorr{i};
end

m_lags = nanmean(lagsBlock,1);
m_corr = nanmean(corrBlock,1);
sem_corr = nanstd(corrBlock,1)./sqrt(i);

%Plot the outputs
fig1 = figure;
% %subplot(2,1,1)
hold on
plot(m_lags,m_corr)
plot(m_lags,m_corr+sem_corr,':')
plot(m_lags,m_corr-sem_corr,':')
title('Correlation Lags Between Song and Neural Activity')
xlabel('Neural Signal Lead (ms)')
ylabel('Correlation')
set(gca,'TickDir','out','Box','off')
xlim([-400,400])

% subplot(2,1,2)
% hold on
% plot(snipLags,m_snipCorr)
% plot(snipLags,m_snipCorr+sem_snipCorr,':')
% plot(snipLags,m_snipCorr-sem_snipCorr,':')
% xlabel('Neural Signal Lead (ms)')
% ylabel('Correlation')
% set(gca,'TickDir','out','Box','off')
% xlim([-100,100])

%Save outputs to file
%save([path 'Over all birds.mat'],'snipLags','m_snipCorr','m_lags','m_corr','peakTime','peakSnip')
% saveas(fig1,[path 'Over all birds plot.fig'])



