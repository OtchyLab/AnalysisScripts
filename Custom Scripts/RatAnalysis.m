function RatAnalysis(logcoords)
%Key variable
numBins = 50;

%Filters
sigma = 1; width = 5;
h = fspecial('gaussian', width, sigma);

%Parse data from pMAT according to Julie's instructions
for i = 1:size(logcoords,1)
%     basePDF(i,:,:) = imfilter(pMAT{i,1},h,'replicate');
%     pbsPDF(i,:,:) = imfilter(pMAT{i,2},h,'replicate');
%     muscPDF(i,:,:) = imfilter(pMAT{i,3},h,'replicate');
%     prePDF(i,:,:) = imfilter(pMAT{i,4},h,'replicate');
    f = figure(100);
    basePDF(i,:,:) = imfilter(ndhist(logcoords{i,1}(:,1), logcoords{i,1}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]),h,'replicate');
    pbsPDF(i,:,:) = imfilter(ndhist(logcoords{i,2}(:,1), logcoords{i,2}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]),h,'replicate');
    muscPDF(i,:,:) = imfilter(ndhist(logcoords{i,3}(:,1), logcoords{i,3}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]),h,'replicate');
    prePDF(i,:,:) = imfilter(ndhist(logcoords{i,4}(:,1), logcoords{i,4}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]),h,'replicate');
    close(f);
%     basePDF(i,:,:) = ndhist(logcoords{i,1}(:,1), logcoords{i,1}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]);
%     pbsPDF(i,:,:) = ndhist(logcoords{i,2}(:,1), logcoords{i,2}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]);
%     muscPDF(i,:,:) = ndhist(logcoords{i,3}(:,1), logcoords{i,3}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]);
%     prePDF(i,:,:) = ndhist(logcoords{i,4}(:,1), logcoords{i,4}(:,2), 'edgesx', linspace(3,8,numBins), 'edgesy', linspace(3,8,numBins), 'prob', 'axis', [3 8 3 8]);
end

%Plot distributons
figure; clf
j = 1;
for i = 1:size(logcoords,1)
    subplot(4,4,j)
    imagesc(squeeze(prePDF(i,:,:))); cs = caxis;
    j = j+1;
    
    subplot(4,4,j)
    imagesc(squeeze(basePDF(i,:,:)), cs);
    j = j+1;
    
    subplot(4,4,j)
    imagesc(squeeze(pbsPDF(i,:,:)), cs);
    j = j+1;
    
    subplot(4,4,j)
    imagesc(squeeze(muscPDF(i,:,:)),cs);
    j = j+1;
end

%Generate the Cost Matrix for EMD computations
COST_MULT_FACTOR= 5;
THRESHOLD = 25;

[R,C] = size(squeeze(basePDF(1,:,:)));
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
%                 D(i,j)= min( [THRESHOLD (COST_MULT_FACTOR*sqrt((10^(r1)-10^(r2))^2+(10^(c1)-10^(c2)^2))] );
            end
        end
    end
end

%Saved cost matrix
% load('CostMatRat-thresh50.mat')

%Calculate all EMDs
for i = 1:size(basePDF,1)
    EMD(i,1) = calcEMD(squeeze(basePDF(i,:,:)), squeeze(prePDF(i,:,:)), D);
    EMD(i,2) = calcEMD(squeeze(basePDF(i,:,:)), squeeze(pbsPDF(i,:,:)), D);
    EMD(i,3) = calcEMD(squeeze(basePDF(i,:,:)), squeeze(muscPDF(i,:,:)), D);
end

%Normalize by pre-base values for shits
for i = 1:size(EMD,2)
    EMDnorm(:,i) = EMD(:,i)./EMD(:,1);
end

%Plot the results
figure(2); clf
subplot(1,2,1); cla
plot([1,2,3;1,2,3;1,2,3;1,2,3], EMD,'ok'); hold on
errorbar([1,2,3], mean(EMD,1), std(EMD,1)./sqrt(size(EMD,1)),'sk', 'LineWidth', 2)
xlim([0.5, 3.5]); %ylim([0,18]);
ylabel('Wasserstein Distance', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10, 'XTick', [1, 2, 3], 'XTickLabels', [{'Control'}, {'PBS'}, {'Musc'}], 'YTick', 0:6:18)

subplot(1,2,2); cla
plot([1,2,3;1,2,3;1,2,3;1,2,3], EMDnorm,'ok'); hold on
errorbar([1,2,3], mean(EMDnorm,1), std(EMDnorm,1)./sqrt(size(EMDnorm,1)),'sk', 'LineWidth', 2)
xlim([0.5, 3.5]); %ylim([0,28]);
ylabel('Normalized Wasserstein Distance', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10, 'XTick', [1, 2, 3], 'XTickLabels', [{'Control'}, {'PBS'}, {'Musc'}], 'YTick', 0:7:28)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4, 3])
%Hypothesis testing
[t_basePBS,p_basePBS] = ttest(EMD(:,1)', EMD(:,2)')
[t_PBSmusc,p_PBSmusc] = ttest(EMD(:,2)', EMD(:,3)')
[t_baseMusc,p_baseMusc] = ttest(EMD(:,1)', EMD(:,3)')

[t_basePBSn,p_basePBSn] = ttest(EMDnorm(:,1)', EMDnorm(:,2)')
[t_PBSmuscn,p_PBSmuscn] = ttest(EMDnorm(:,2)', EMDnorm(:,3)')
[t_baseMuscn,p_baseMuscn] = ttest(EMDnorm(:,1)', EMDnorm(:,3)')

% p_basePBS = ranksum(EMD(:,1)', EMD(:,2)')
% p_baseMusc = ranksum(EMD(:,1)', EMD(:,3)')
% 
% p_basePBSn = ranksum(EMDnorm(:,1)', EMDnorm(:,2)')
% p_baseMuscn = ranksum(EMDnorm(:,1)', EMDnorm(:,3)')


a = 1;
