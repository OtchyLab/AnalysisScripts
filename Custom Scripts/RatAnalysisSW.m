function RatAnalysisSW(logcoords)
%Processes the log-log data from Steffen's opto experiment. Input data ('logcoords') is a cell array of size (mxn):
%       m = individual animals (rows)
%       n = experimental treatments (columns)
%               the first column is the pooled 'light on' data
%               all subsequent columns are the individual sessions that are 'light off''
%
%The general steps for this function are:
%       1) generate 2-D pdfs of log-log Data, apply simple gaussian filter, and plot
%       2) Calculate the cost/distance matrix for the pdfs (to be used for EMD measurements)
%       3) Calculate the EMD/Wasserstein distance between all 'lights off' conditions
%       4) Calculate the EMD/Wasserstein distance between all 'lights off' conditions and the 'lights on' condition
%       5) Plot summary figures of EMDs
%       6) Hypothesis testing via paired t-test

%Key variable
numBins = 50;

%Filters
sigma = 1; width = 5;
h = fspecial('gaussian', width, sigma);

%Parse data from pMAT according to Steffen's instructions
for i = 1:size(logcoords,1)
    %temp figure to avoid fucking up real figures
    f = figure;
    onPDF(i,:,:) = imfilter(ndhist(logcoords{i,1}(:,1), logcoords{i,1}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8])',h,'replicate');
    
    for j = 2:size(logcoords,2)
        offPDF(i,j-1,:,:) = imfilter(ndhist(logcoords{i,j}(:,1), logcoords{i,j}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8])',h,'replicate');
    end
    close(f);
end

%Plot distributons
for i = 1:size(logcoords,1)
    subs = size(logcoords,2);
    h(i) = figure(i); clf
    subplot(1,subs,1)
    imagesc(squeeze(onPDF(i,:,:))'); cs = caxis; axis xy
    xlabel('log(Tap1 Dur)'); ylabel('log(IPI-Tap1 Dur)')
    set(gca,'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'XTick', linspace(1,numBins,3), 'YTick', linspace(1,numBins,3),'XTickLabel', 3:2:8, 'YTickLabel',3:2:8)
    
    for j = 2:subs;
        subplot(1,subs,j)
        imagesc(squeeze(offPDF(i,j-1,:,:)), cs); axis xy
        set(gca,'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'XTick', linspace(1,numBins,3), 'YTick', linspace(1,numBins,3),'XTickLabel', 3:2:8, 'YTickLabel',3:2:8)
    end
    set(gcf,'Units','Inches','Position',[0,0,10,1.5])
end

%Generate the Cost Matrix for EMD computations
COST_MULT_FACTOR= 5;
THRESHOLD = 25;

[R,C] = size(squeeze(onPDF(1,:,:)));
D= zeros(R*C,R*C,'int32');
j= 0;
for c1=1:C
    for r1=1:R
        j= j+1;
        i= 0;
        for c2=1:C
            for r2=1:R
                i = i+1;
                D(i,j)= min( [THRESHOLD (COST_MULT_FACTOR*sqrt((r1-r2)^2+(c1-c2)^2))] );
            end
        end
    end
end

%Saved cost matrix
% load('CostMatRat-thresh50.mat')

%Calculate all EMD between lights on and all lights off
onoffEMD = [];
for i = 1:size(onPDF,1)
    for j = 1:size(squeeze(offPDF(i,:,:,:)),1)
        onoffEMD(i,j) = calcEMD(squeeze(onPDF(i,:,:)), squeeze(offPDF(i,j,:,:)), D);
    end
end

%Calculate all EMD between all lights off
offoffEMD = NaN(size(onPDF,1),size(squeeze(offPDF(i,:,:,:)),1),size(squeeze(offPDF(i,:,:,:)),1));
for i = 1:size(onPDF,1)
    for j = 1:size(squeeze(offPDF(i,:,:,:)),1)
        for k = 1:size(squeeze(offPDF(i,:,:,:)),1)
            if j ~= k
                offoffEMD(i,j,k) = calcEMD(squeeze(offPDF(i,j,:,:)), squeeze(offPDF(i,k,:,:)), D);
            end
        end
    end
end

%Means for each animal
onoffEMDm = mean(onoffEMD,2);
offoffEMDm = mean(nanmean(offoffEMD,3),2);

%Plot summary figure
h(size(logcoords,1)+1) = figure; clf
plot(1*ones(1,length(offoffEMDm)), offoffEMDm, '.k','MarkerSize', 10); hold on
plot(2*ones(1,length(onoffEMDm)), onoffEMDm, '.r','MarkerSize', 10);
errorbar([1,2],[mean(offoffEMDm), mean(onoffEMDm)],[std(offoffEMDm), std(onoffEMDm)],'sk', 'LineWidth',1.25)
xlim([0.5, 2.5]); ylim([4,14])
ylabel('Wasserstein Distance')
set(gca,'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'XTick', [1,2],'XTickLabel', [{'Control'},{'Light ON'}], 'YTick', 4:4:16,'YTickLabel', 4:4:16)
set(gcf,'Units','Inches','Position',[0,0,4,3])

%Hypothesis Testing
[t_, p_] = ttest(onoffEMDm,offoffEMDm)


save('C:\Users\Tim\Desktop\Nif Project Figures\Steffen Analysis\First Pass Analysis.mat')
