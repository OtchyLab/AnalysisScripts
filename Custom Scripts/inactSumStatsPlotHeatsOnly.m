function inactSumStatsPlotHeatsOnly(inactStats)

summary = [];

name = getStructField(inactStats,'name');
date = getStructField(inactStats,'date');
type = getStructField(inactStats,'type');
% d2emd = getStructField(inactStats,'D2emd');

preY = getStructField(inactStats,'preY');
injY = getStructField(inactStats,'injY');
postY = getStructField(inactStats,'postY');

%use the smoothed heatmaps
preN = getStructField(inactStats,'preNsm');
injN = getStructField(inactStats,'injNsm');
postN = getStructField(inactStats,'postNsm');

birds = unique(name);
conditions = [{'inact'}, {'PBS'}, {'elev'}];
color = {'r', 'b', 'g'};

%Assemble figures for individuall birds
for i = 1:length(birds)
    birdMask = strcmp(birds(i), name);
    
    %Bar plot all of the similarity/distance measures on a single figure
    h(3*i-2) = figure(3*i-2); clf
    
%     for j = 1:length(conditions)
%         typeMask = strcmp(conditions(j), type);
%         
%         summary(end+1).name = birds{i};
%         summary(end).type = char(conditions(j));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% 2-D measures %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         %Bar plot the 2D dist
%         subSet = d2emd((birdMask & typeMask),:);
%         subplot(2,7,9)
%         hold on
%         for k = 1:length(subSet(:,3))
%             plot([1,j+1], [subSet(k,3), subSet(k,1)], 'DisplayName', char(conditions(j)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o');
%         end
%                     %format
%         xlim([0.5,4.5]); %ylim([0,1.1]); %axis square
%         ylabel('Wasserstein Distance', 'FontSize', 8)
%         title([char(birds(i)) ' D2 Wasserstein Dist'], 'FontSize', 8);
%         set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'LineWidth', 2, 'FontSize', 10)
%         
%         summary(end).D2emdmean = mean([subSet(:,3), subSet(:,1)],1);
%         summary(end).D2emdstd = std([subSet(:,3), subSet(:,1)],1);
% 
%         %Do the day-over-day euclidean dist measure
%         subDate = squeeze(date(birdMask));
%         preNs = squeeze(preN(birdMask, :, :));
%         pd = [];
%         for m = 1:size(preNs,1)
%             for n = 1:size(preNs,1)
%                 d1 = str2num(subDate{m}((end-1):end));
%                 d2 = str2num(subDate{n}((end-1):end));
%                 if (d2-d1) == 1 && d1 ~=2 && d2 ~=2 %Fixed bad file
%                     t1 = squeeze(preNs(m,:,:));
%                     t2 = squeeze(preNs(n,:,:));
%                     pd = [pd; calcEMD(t1,t2)]; %Stack
%                 end
%             end
%         end
%         pd = pd(pd~=0);
%         summary(end).D2distDay = nanmean(pd);
% 
%     end
%     %Set figure size
%     set(gcf, 'Units', 'Inches', 'Position', [0 0 15 4])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  Composite Distribution Images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %Pre & post
    preYmean = mean(preY(birdMask, :),1);
    postYmean = mean(postY(birdMask, :),1);
    
    preNmean = squeeze(mean(preN(birdMask, :, :),1));
    postNmean = squeeze(mean(postN(birdMask, :, :),1));
    
    %Generate unique masks for conditions
    inactMask = strcmp(conditions(1), type);
    PBSMask = strcmp(conditions(2), type);
    elevMask = strcmp(conditions(3), type);
    
    %Treatment means
    inactYmean = mean(injY(birdMask & inactMask, :),1);
    PBSYmean = mean(injY(birdMask & PBSMask, :),1);
    elevYmean = mean(injY(birdMask & elevMask, :),1);
    
    inactNmean = squeeze(mean(injN(birdMask & inactMask, :, :),1));
    PBSNmean = squeeze(mean(injN(birdMask & PBSMask, :, :),1));
    elevNmean = squeeze(mean(injN(birdMask & elevMask, :, :),1));
    
    h(3*i-1) = figure(3*i-1); clf
    plot(0:1.25:350, preYmean, 'LineWidth', 2, 'DisplayName', 'Pre')
    hold on
    plot(0:1.25:350, inactYmean, 'LineWidth', 2, 'DisplayName', 'Inact')
    plot(0:1.25:350, PBSYmean, 'LineWidth', 2, 'DisplayName', 'PBS')
    plot(0:1.25:350, elevYmean, 'LineWidth', 2, 'DisplayName', 'Elev')
    plot(0:1.25:350, postYmean, 'LineWidth', 2, 'DisplayName', 'Post')

    %Format the figure
    axis tight
    xlim([0 400])
    ylim([0,max([preYmean, postYmean])*1.1])
    box off
    set(gca, 'TickDir', 'out')
    xlabel('Syl Dur (ms)', 'FontSize', 16)
    ylabel('P(t)', 'FontSize', 16)
    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 5 5])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    title([char(birds(i)) ' D1 PDF'], 'FontSize', 16);
    legend('boxoff')
    
    h(3*i) = figure(3*i); clf;  hold on
    subplot(5,1,1)
    imagesc(preNmean'); sc = caxis; colormap(jet); axis xy
    box off
    set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    title([char(birds(i)) ' D2 PDF'], 'FontSize', 16);
    
    subplot(5,1,2)
    imagesc(inactNmean', sc); colormap(jet); axis xy
    box off
    set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    
     subplot(5,1,3)
    imagesc(PBSNmean', sc); colormap(jet); axis xy
    box off
    set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    
    subplot(5,1,4)
    imagesc(elevNmean', sc); colormap(jet); axis xy
    box off
    set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [], 'XTick', [1:20:79], 'XTickLabels', [])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    
    subplot(5,1,5)
    imagesc(postNmean', sc); colormap(jet); axis xy
    box off
    set(gca, 'TickDir', 'out', 'YTick', [1:20:41], 'YTickLabels', [-4, -2, 0], 'XTick', [1:20:79], 'XTickLabels', [0, 100, 200, 300])
    set(gca, 'LineWidth', 3, 'FontSize', 16)
    xlabel('Syl Dur (ms)'); ylabel('log(Entropy)')

    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 4 10])

end
 
%Extract the data from the summary structure
Sumname = getStructField(summary,'name');
Sumtype = getStructField(summary,'type');

SumD2emdmean = getStructField(summary,'D2emdmean');
SumD2distDay = unique(getStructField(summary,'D2distDay'));

%Calculate summary stats and plot results
sumbirds = unique(Sumname);
sumconditions = [{'inact'}, {'PBS'}, {'elev'}];

% h(end+1) = figure(20); clf
% for i = 1:length(sumconditions)
%     %Set the condition mask
%     typeMask = strcmp(sumconditions(i), Sumtype); 
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%% 2-D measures %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Scatter plot D2 Cosine with errorbars
%     subSet = SumD2emdmean(typeMask,:);
%     
%     subplot(2,7,10)
%     hold on
%     for j = 1:length(subSet(:,1))
%         plot([1,i+1], [subSet(j,1), subSet(j,2)], 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);  
%     end
%     errorbar(i+1, nanmean(subSet(:,2)), nanstd(subSet(:,2)), 'Color', color{i}, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);    
%         D2emdPP(:,i) = subSet(:,1);
%     if i ==3
%         errorbar(1, nanmean(D2emdPP(:)), nanstd(D2emdPP(:)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);
%     end
%         %format
%     xlim([0.5,4.5]); %ylim([0,1.1]); %axis square
%     ylabel('Wasserstein Distance', 'FontSize', 10)
%     title(['D2 EMD'], 'FontSize', 8);
%     set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'LineWidth', 2, 'FontSize', 10)    
%     
% end
% set(gcf, 'Units', 'Inches', 'Position', [0 0 15 4])
% 
% %Generate
% ppDist = [];
% inactDist = [];
% PBSDist = [];
% elevDist = [];
% 
% %Organize data for plotting
% ddDist = SumD2distDay;
% for j = 1:length(sumbirds)
%     birdMask = strcmp(sumbirds(j), Sumname);
%     
%     subdist = SumD2emdmean(birdMask,:);
%     
%     ppDist = [ppDist; mean(subdist(:,1))];
%     
%     for i = 1:length(sumconditions)
%         %Set the condition mask
%         typeMask = strcmp(sumconditions(i), Sumtype);
% 
%         %Grab subsets for this treatment
%         subdist = SumD2emdmean(birdMask & typeMask,:);
% 
%         %Parse
%         if i == 1
%             inactDist = [inactDist; subdist(:,2)];
%         elseif i == 2
%             PBSDist = [PBSDist; subdist(:,2)];        
%         elseif i ==3
%             elevDist = [elevDist; subdist(:,2)];        
%         end
% 
%     end
% end
% 
% %Normalize to baseline variation/difference
% % for i = 1:size(ddKL,1)
% %     ppKL(i) = ppKL(i)/ddKL(i);
% %     inactKL(i) = inactKL(i)/ddKL(i);
% %     PBSKL(i) = PBSKL(i)/ddKL(i);
% %     elevKL(i) = elevKL(i)/ddKL(i);
% %     ddKL(i) = ddKL(i)/ddKL(i);
% %     
% %     ppDist(i) = ppDist(i)/ddDist(i);
% %     inactDist(i) = inactDist(i)/ddDist(i);
% %     PBSDist(i) = PBSDist(i)/ddDist(i);
% %     elevDist(i) = elevDist(i)/ddDist(i);
% %     ddDist(i) = ddDist(i)/ddDist(i);
% % end
% 
% % close all
% 
% %Plot polar plot of Summary differences in 2 metrics
% h(end+1) = figure(101); clf
% % subplot(1,2,1)
% % p(1) = plot(1*ones(length(ddKL),1), ddKL, '.k'); hold on
% % p(2) = plot(2*ones(length(ppKL),1), ppKL, '.k'); 
% % p(3) = plot(3*ones(length(inactKL),1), inactKL, '.k');
% % p(4) = plot(4*ones(length(PBSKL),1), PBSKL, '.k');
% % p(5) = plot(5*ones(length(elevKL),1), elevKL, '.k');
% % eb(1) = errorbar(1:5, [mean(ddKL), mean(ppKL), mean(inactKL), mean(PBSKL), mean(elevKL)], [std(ddKL), std(ppKL), std(inactKL), std(PBSKL), std(elevKL)]./sqrt(5),'sk');
% % xlim([0.5,5.5]); ylim([0.75,4.5]); %axis square
% % ylabel('Norm Divergence (AU)', 'FontSize', 12)
% % title(['KL Divergence'], 'FontSize', 12);
% % set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', [{'Cntl'}, {'P-P'}, conditions], 'YTick', [1,2:4.5], 'LineWidth', 2, 'FontSize', 10)
%     
% % subplot(1,2,2)
% q(1) = plot(1*ones(length(ddDist),1), ddDist, '.k'); hold on
% q(2) = plot(2*ones(length(ppDist),1), ppDist, '.k');
% q(3) = plot(3*ones(length(inactDist),1), inactDist, '.k');
% q(4) = plot(4*ones(length(PBSDist),1), PBSDist, '.k');
% q(5) = plot(5*ones(length(elevDist),1), elevDist, '.k');
% eb(2) = errorbar(1:5, [mean(ddDist), mean(ppDist), mean(inactDist), mean(PBSDist), mean(elevDist)], [std(ddDist), std(ppDist), std(inactDist), std(PBSDist), std(elevDist)]./sqrt(5),'sk');
% xlim([0.5,5.5]); ylim([0,25]); %axis square
% ylabel('Norm EMD (AU)', 'FontSize', 12)
% title(['Earth Mover Distance'], 'FontSize', 12);
% set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', [{'Cntl'}, {'P-P'}, conditions], 'YTick', [1:4:8], 'LineWidth', 2, 'FontSize', 10)   
% 
% for i = 1:5
%     subplot(1,2,1)
%     set(p(i), 'Marker', '.', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'LineStyle', 'none');
%     set(eb(1), 'Marker', 's', 'Color', 'k', 'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5);
%     
%     subplot(1,2,2)
%     set(q(i), 'Marker', '.', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'LineStyle', 'none');
%     set(eb(2), 'Marker', 's', 'Color', 'k', 'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5);
% end
% 
% set(gcf, 'Units', 'Inches', 'Position', [0 0 6 3])
% 
% %Hypothesis Testing
% [testOut.ddpp,pOut.ddpp,ci.ddpp,~] = ttest([ddDist],[ppDist]);
% [testOut.ddInact,pOut.ddInact,ci.ddInact,~] = ttest([ddDist],[inactDist]);
% [testOut.ddPBS,pOut.ddPBS,ci.ddPBS,~] = ttest([ddDist],[PBSDist]);
% [testOut.ddElev,pOut.ddElev,ci.ddElev,~] = ttest([ddDist],[elevDist]);

%Save figure and data
saveFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\';
saveName = 'Summary Inactivation Data 05152015 temp';

%Save figures
% savefig(h, [saveFolder, saveName '.fig']);
% close all; clear('h')

%Save summary data
save([saveFolder, saveName '.mat']);






