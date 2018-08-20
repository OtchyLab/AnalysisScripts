function inactSumStatsPlot(inactStats)

summary = [];

name = getStructField(inactStats,'name');
type = getStructField(inactStats,'type');
d1corr = getStructField(inactStats,'D1corr');
d1dist = getStructField(inactStats,'D1dist');
d2corr = getStructField(inactStats,'D2corr');
d2dist = getStructField(inactStats,'D2dist');

birds = unique(name);
% conditions = unique(type);
conditions = [{'Inact'}, {'PBS'}, {'Elev'}];
m = 1;
color = {'r', 'b', 'g'};
for i = 1:length(birds)
    birdMask = strcmp(birds(i), name);
    
    for j = 1:length(conditions)
        typeMask = strcmp(conditions(j), type);
        
        subSet = d1corr((birdMask & typeMask),:);
        %Scatter plot the 1D correlations
        h(2*i-1) = figure(2*i-1);
        subplot(2,2,1)
        hold on
        scatter(subSet(:,3), subSet(:,1), 'DisplayName', char(conditions(j)), 'MarkerEdgeColor', color{j});
        
        %Bar plot the 1D correlations
        h(2*i) = figure(2*i);
        subplot(2,2,1)
        hold on
        for k = 1:length(subSet(:,3))
            plot([1,j+1], [subSet(k,3)./subSet(k,3), subSet(k,1)./subSet(k,3)], 'DisplayName', char(conditions(j)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o');
        end
        
        summary(end+1).name = birds(i);
        summary(end).type = char(conditions(j));
        summary(end).labels = [{'pre-post'}; {'pre-inj'}];
        summary(end).D1corrmean = mean([subSet(:,3), subSet(:,1)],1);
        summary(end).D1corrstd = std([subSet(:,3), subSet(:,1)],1);
        
        subSet = d1dist((birdMask & typeMask),:);
        %Scatter plot the 1D correlations
        h(2*i-1) = figure(2*i-1);
        subplot(2,2,2)
        hold on
        scatter(subSet(:,3), subSet(:,1), 'DisplayName', char(conditions(j)), 'MarkerEdgeColor', color{j});

        %Bar plot the 1D dist
        h(2*i) = figure(2*i);
        subplot(2,2,2)
        hold on
        for k = 1:length(subSet(:,3))
            plot([1,j+1], [subSet(k,3)./subSet(k,3), subSet(k,1)./subSet(k,3)], 'DisplayName', char(conditions(j)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o');
        end
        
        summary(end).D1distmean = mean([subSet(:,3), subSet(:,1)],1);
        summary(end).D1diststd = std([subSet(:,3), subSet(:,1)],1);
        
        subSet = d2corr((birdMask & typeMask),:);
        %Scatter plot the 1D correlations
        h(2*i-1) = figure(2*i-1);
        subplot(2,2,3)
        hold on
        scatter(subSet(:,3), subSet(:,1), 'DisplayName', char(conditions(j)), 'MarkerEdgeColor', color{j});
        
        %Bar plot the 2D correlations
        h(2*i) = figure(2*i);
        subplot(2,2,3)
        hold on
        for k = 1:length(subSet(:,3))
            plot([1,j+1], [subSet(k,3)./subSet(k,3), subSet(k,1)./subSet(k,3)], 'DisplayName', char(conditions(j)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o');
        end        
        
        summary(end).D2corrmean = mean([subSet(:,3), subSet(:,1)],1);
        summary(end).D2corrstd = std([subSet(:,3), subSet(:,1)],1);
        
        subSet = d2dist((birdMask & typeMask),:);
        %Scatter plot the 1D correlations
        h(2*i-1) = figure(2*i-1);
        subplot(2,2,4)
        hold on
        scatter(subSet(:,3), subSet(:,1), 'DisplayName', char(conditions(j)), 'MarkerEdgeColor', color{j});
        
        %Bar plot the 2D dist
        h(2*i) = figure(2*i);
        subplot(2,2,4)
        hold on
        for k = 1:length(subSet(:,3))
            plot([1,j+1], [subSet(k,3)./subSet(k,3), subSet(k,1)./subSet(k,3)], 'DisplayName', char(conditions(j)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o');
        end        

        summary(end).D2distmean = mean([subSet(:,3), subSet(:,1)],1);
        summary(end).D2diststd = std([subSet(:,3), subSet(:,1)],1);
        
    end
    
    h(2*i-1) = figure(2*i-1);
    subplot(2,2,1)
    plot([0,1], [0,1], '--k')
    xlim([0,1]); ylim([0,1]); axis square
    xlabel('Pre-Post Correlation', 'FontSize', 10)
    ylabel('Pre-Inj Correlation', 'FontSize', 10)
    title([char(birds(i)) ' D1 Correlation'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.5, 1], 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,2)
    plot([0,1], [0,1], '--k')
    xlim([0,0.25]); ylim([0,0.25]); axis square
    xlabel('Pre-Post Distance', 'FontSize', 10)
    ylabel('Pre-Inj Distance', 'FontSize', 10)
    title([char(birds(i)) ' D1 Distance'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.1, 0.2], 'YTick', [0, 0.1, 0.2], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,3)
    plot([0,1], [0,1], '--k')
    xlim([0,1]); ylim([0,1]); axis square
    xlabel('Pre-Post Correlation', 'FontSize', 10)
    ylabel('Pre-Inj Correlation', 'FontSize', 10)
    title([char(birds(i)) ' D2 Correlation'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.5, 1], 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,4)
    plot([0,1], [0,1], '--k')
    xlim([0,0.05]); ylim([0,0.05]); axis square
    xlabel('Pre-Post Distance', 'FontSize', 10)
    ylabel('Pre-Inj Distance', 'FontSize', 10)
    title([char(birds(i)) ' D2 Distance'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.2, 0.4], 'YTick', [0, 0.02, 0.04], 'LineWidth', 2, 'FontSize', 10)
    
    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 5 5])

    legend(char(conditions));
    
    h(2*i) = figure(2*i);
    subplot(2,2,1)
    plot([0,5], [1,1], ':k')
    xlim([0.5,4.5]); ylim([0,1.5]); axis square
    ylabel('Norm. Correlation', 'FontSize', 10)
    title([char(birds(i)) ' D1 Correlation'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 0.5, 1, 1.5], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,2)
    plot([0,5], [1,1], ':k')
    xlim([0.5,4.5]); ylim([0,10]); axis square
    ylabel('Norm. Distance', 'FontSize', 10)
    title([char(birds(i)) ' D1 Distance'], 'FontSize', 10);
        set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 2.5, 5, 7.5, 10], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,3)
    plot([0,5], [1,1], ':k')
    xlim([0.5,4.5]); ylim([0,1.5]); axis square
    ylabel('Norm. Correlation', 'FontSize', 10)
    title([char(birds(i)) ' D2 Correlation'], 'FontSize', 10);
        set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 0.5, 1, 1.5], 'LineWidth', 2, 'FontSize', 10)

    subplot(2,2,4)
    plot([0,5], [1,1], ':k')
    xlim([0.5,4.5]); ylim([0,5]); axis square
    ylabel('Norm. Distance', 'FontSize', 10)
    title([char(birds(i)) ' D2 Distance'], 'FontSize', 10);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 2.5, 5], 'LineWidth', 2, 'FontSize', 10)
    
    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 5 5])

%     legend(char(conditions)); legend('boxoff')
end

        
Sumname = getStructField(summary,'name');
Sumtype = getStructField(summary,'type');
SumD1corrmean = getStructField(summary,'D1corrmean');
SumD1corrstd = getStructField(summary,'D1corrstd');
SumD1distmean = getStructField(summary,'D1distmean');
SumD1diststd = getStructField(summary,'D1diststd');
SumD2corrmean = getStructField(summary,'D2corrmean');
SumD2corrstd = getStructField(summary,'D2corrstd');
SumD2distmean = getStructField(summary,'D2distmean');
SumD2diststd = getStructField(summary,'D2diststd');

sumbirds = unique(Sumname);
sumconditions = [{'Inact'}, {'PBS'}, {'Elev'}];

for i = 1:length(sumconditions)
    typeMask = strcmp(sumconditions(i), Sumtype);
    
    %D1 Correlations
    subSet = SumD1corrmean(typeMask,:);
    subStd = SumD1corrstd(typeMask,:);
    
    figure(100);
    subplot(2,2,1)
    hold on
    errorbarxy(subSet(:,1), subSet(:,2), subStd(:,1), subStd(:,2),{['o' color{i}],color{i},color{i}});

    figure(101)
    subplot(2,2,1)
    hold on
    for j = 1:length(subSet(:,1))
        plot([1,i+1], [subSet(j,1)./subSet(j,1), subSet(j,2)./subSet(j,1)], 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);  
    end
    errorbar(i+1, nanmean(subSet(:,2)./subSet(:,1)), nanstd(subSet(:,2)./subSet(:,1)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);
    
    %D1 Distance
    subSet = SumD1distmean(typeMask,:);
    subStd = SumD1diststd(typeMask,:);
    
    figure(100);
    subplot(2,2,2)
    hold on
    errorbarxy(subSet(:,1), subSet(:,2), subStd(:,1), subStd(:,2),{['o' color{i}],color{i},color{i}});

    figure(101)
    subplot(2,2,2)
    hold on
    for j = 1:length(subSet(:,1))
        plot([1,i+1], [subSet(j,1)./subSet(j,1), subSet(j,2)./subSet(j,1)], 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);  
    end
    errorbar(i+1, nanmean(subSet(:,2)./subSet(:,1)), nanstd(subSet(:,2)./subSet(:,1)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);    
    
    %D2 Correlations
    subSet = SumD2corrmean(typeMask,:);
    subStd = SumD2corrstd(typeMask,:);
    
    figure(100);
    subplot(2,2,3)
    hold on
    errorbarxy(subSet(:,1), subSet(:,2), subStd(:,1), subStd(:,2),{['o' color{i}],color{i},color{i}});

    figure(101)
    subplot(2,2,3)
    hold on
    for j = 1:length(subSet(:,1))
        plot([1,i+1], [subSet(j,1)./subSet(j,1), subSet(j,2)./subSet(j,1)], 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);  
    end
    errorbar(i+1, nanmean(subSet(:,2)./subSet(:,1)), nanstd(subSet(:,2)./subSet(:,1)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);    
    
    %D2 Distance
    subSet = SumD2distmean(typeMask,:);
    subStd = SumD2diststd(typeMask,:);
    
    figure(100);
    subplot(2,2,4)
    hold on
    errorbarxy(subSet(:,1), subSet(:,2), subStd(:,1), subStd(:,2),{['o' color{i}],color{i},color{i}});

    figure(101)
    subplot(2,2,4)
    hold on
    for j = 1:length(subSet(:,1))
        plot([1,i+1], [subSet(j,1)./subSet(j,1), subSet(j,2)./subSet(j,1)], 'Color', [.5, .5, .5], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);  
    end
    errorbar(i+1, nanmean(subSet(:,2)./subSet(:,1)), nanstd(subSet(:,2)./subSet(:,1)), 'Color', 'k', 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 8);      
    
end

h(end+1) = figure(100);
subplot(2,2,1)
hold on
plot([0,1], [0,1], '--k')
xlim([0,1]); ylim([0,1]); axis square
xlabel('Pre-Post Correlation', 'FontSize', 10)
ylabel('Pre-Inj Correlation', 'FontSize', 10)
title(['Summary D1 Correlation'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.5, 1], 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,2)
hold on
plot([0,2], [0,2], '--k')
xlim([0,0.4]); ylim([0,0.4]); axis square
xlabel('Pre-Post Distance', 'FontSize', 10)
ylabel('Pre-Inj Distance', 'FontSize', 10)
title(['Summary D1 Distance'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.2, 0.4], 'YTick', [0, 0.2, 0.4], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,3)
hold on
plot([0,1], [0,1], '--k')
xlim([0,1]); ylim([0,1]); axis square
xlabel('Pre-Post Correlation', 'FontSize', 10)
ylabel('Pre-Inj Correlation', 'FontSize', 10)
title(['Summary D2 Correlation'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.5, 1], 'YTick', [0, 0.5, 1], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,4)
hold on
plot([0,1], [0,1], '--k')
xlim([0,0.05]); ylim([0,0.05]); axis square
xlabel('Pre-Post Distance', 'FontSize', 10)
ylabel('Pre-Inj Distance', 'FontSize', 10)
title(['Summary D2 Distance'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 0.2, 0.4], 'YTick', [0, 0.02, 0.04], 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 5])

% legend(char(conditions));

h(end+1) = figure(101);
subplot(2,2,1)
hold on
plot([0,5], [1,1], ':k')
xlim([0.5,4.5]); ylim([0,1.5]); axis square
ylabel('Norm. Correlation', 'FontSize', 10)
title(['Summary D1 Correlation'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 0.5, 1, 1.5], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,2)
hold on
plot([0,5], [1,1], ':k')
xlim([0.5,4.5]); ylim([0,10]); axis square
ylabel('Norm. Distance', 'FontSize', 10)
title(['Summary D1 Distance'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 1, 5, 10], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,3)
hold on
plot([0,5], [1,1], ':k')
xlim([0.5,4.5]); ylim([0,1.5]); axis square
ylabel('Norm. Correlation', 'FontSize', 10)
title(['Summary D2 Correlation'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 0.5, 1, 1.5], 'LineWidth', 2, 'FontSize', 10)

subplot(2,2,4)
hold on
plot([0,5], [1,1], ':k')
xlim([0.5,4.5]); ylim([0,3]); axis square
ylabel('Norm. Distance', 'FontSize', 10)
title(['Summary D2 Distance'], 'FontSize', 10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'P-P'}, conditions], 'YTick', [0, 1, 2, 3], 'LineWidth', 2, 'FontSize', 10)

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 5])

%Save figure and data
saveFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\';
saveName = 'Summary Inactivation Data';

%Save summary data
save([saveFolder, saveName '.mat'], 'summary');

%Save figures
savefig(h, [saveFolder, saveName '.fig']);


