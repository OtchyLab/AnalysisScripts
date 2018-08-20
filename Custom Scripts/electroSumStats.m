function electroSumStats(longStats)
%Takes the output of LongSong.m
summary = [];

name = getStructField(longStats,'name');
type = getStructField(longStats,'type');
durCorr = getStructField(longStats,'durCorr');
xs = getStructField(longStats,'xs');
SCsnipsM = getStructField(longStats,'SCsnipsM');
SCsnipsS = getStructField(longStats,'SCsnipsS');
MFsnipsM = getStructField(longStats,'MFsnipsM');
MFsnipsS = getStructField(longStats,'MFsnipsS');
MCsnipsM = getStructField(longStats,'MCsnipsM');
MCsnipsS = getStructField(longStats,'MCsnipsS');

birds = unique(name);
conditions = [{'hit'}];

color = {'r', 'b', 'g'};
for i = 1:length(conditions)
    typeMask = strcmp(conditions(i), type);

    %Scatterplot the recovery of corelation
    subSet = durCorr((typeMask),:);
    x = xs(typeMask);
    
    h(1) = figure(1); clf
    subplot(2,2,1:2); cla
    hold on
    %Real data
    for j = 1:size(subSet,1)
        scatter(x(j).D, subSet(j, x(j).D), 'sb')
        scatter(x(j).N, subSet(j, x(j).N), 'sk')
        scatter(x(j).L, subSet(j, x(j).L), 'sr')        
    end
    errorbar(x(j).D, mean(subSet(:,x(j).D)), std(subSet(:,x(j).D)), 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(x(j).N, mean(subSet(:,x(j).N)), std(subSet(:,x(j).N)), 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    errorbar(x(j).L, mean(subSet(:,x(j).L)), std(subSet(:,x(j).L)), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    plot(x(j).all, mean(subSet(:, x(j).all)), '--k')
    
    %Fake control data
    fakex = [1, 8, 9;1, 8, 9;1, 8, 9];
    fakey = [1, 0.92, 0.93; 1, 0.93, 0.91; 1, 0.88, 0.91];
    for j = 1:size(subSet,1)
        scatter(fakex(j, 1:2), fakey(j, 1:2), '.b')
        scatter(fakex(j, 3), fakey(j, 3), '.k')      
    end
    errorbar(fakex(j, 1:2), mean(fakey(:, 1:2),1), std(fakey(:, 1:2),1), 'db', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(fakex(j, 3), mean(fakey(:,3), 1), std(fakey(:,3), 1), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    plot(fakex(j, :), mean(fakey, 1), ':k')
    xlim([0.5, 9.5])
    ylim([0.4,1.1]);
    ylabel('Recovery of Duration', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', x(j).all, 'XTickLabels', [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}])
    
    subplot(2,2,3); cla
    hold on
    day = 2:2:8; night = 3:2:7;
    recovDay = subSet(:, day+1) - subSet(:, day);
    recovNight = subSet(:, night+1) - subSet(:, night);
    scatter(ones(length(recovDay(:)),1), recovDay(:), '.b')
    scatter(2*ones(length(recovNight(:)),1), recovNight(:), '.k')
    errorbar([1], [mean(recovDay(:))], [std(recovDay(:))], 'ob', 'MarkerSize', 8)
    errorbar([2], [mean(recovNight(:))], [std(recovNight(:))], 'ok', 'MarkerSize', 8)
    plot([1, 2], [mean(recovDay(:)), mean(recovNight(:))], '--k', 'LineWidth', 1)
    plot([0.5, 2.5], [0,0], ':k')
    xlim([0.5, 2.5])
    ylim([-0.3,0.3]);
    ylabel('Change in Corr', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', [1, 2], 'XTickLabels', [{'Over Day'}, {'Over Night'}], 'YTick', [-.3, 0, .3])
    
%     summary(end+1).name = birds;
%     summary(end).type = char(type(j));
%     summary(end). = [{'pre-post'}; {'pre-inj'}];
%     summary(end).D1corrmean = mean([subSet(:,3), subSet(:,1)],1);
%     summary(end).D1corrstd = std([subSet(:,3), subSet(:,1)],1);

%Scatterplot the recovery of Motif Fraction
    subSet = MFsnipsM((typeMask),:);
    x = xs(typeMask);
    
    h(2) = figure(2); clf
    subplot(2,2,1:2); cla
    hold on
    %Real data
    for j = 1:size(subSet,1)
        scatter(x(j).D, subSet(j, x(j).D)/subSet(j, 1), 'sb')
        scatter(x(j).N, subSet(j, x(j).N)/subSet(j, 1), 'sk')
        scatter(x(j).L, subSet(j, x(j).L)/subSet(j, 1), 'sr')      
        t(j,:) = subSet(j,:)./subSet(j, 1);
    end
    errorbar(x(j).D, mean(t(:,x(j).D)), std(t(:,x(j).D)), 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(x(j).N, mean(t(:,x(j).N)), std(t(:,x(j).N)), 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    errorbar(x(j).L, mean(t(:,x(j).L)), std(t(:,x(j).L)), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    plot(x(j).all, mean(t(:, x(j).all)), '--k')
    
    %Fake control data
    fakex = [1, 8, 9;1, 8, 9;1, 8, 9];
    fakey = [1, 0.92, 0.93; 1, 1.02, 1; 1, 0.93, 0.91];
    for j = 1:size(subSet,1)
        scatter(fakex(j, 1:2), fakey(j, 1:2), '.b')
        scatter(fakex(j, 3), fakey(j, 3), '.k')      
    end
    errorbar(fakex(j, 1:2), mean(fakey(:, 1:2),1), std(fakey(:, 1:2),1), 'db', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(fakex(j, 3), mean(fakey(:,3), 1), std(fakey(:,3), 1), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    plot(fakex(j, :), mean(fakey, 1), ':k')
    xlim([0.5, 9.5])
    ylim([0.2,1.1]);
    ylabel('Recovery of Motif Frac', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', x(j).all, 'XTickLabels', [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}])
    
    subplot(2,2,3); cla
    hold on
    day = 2:2:8; night = 3:2:7;
    recovDay = subSet(:, day+1) - subSet(:, day);
    recovNight = subSet(:, night+1) - subSet(:, night);
    scatter(ones(length(recovDay(:)),1), recovDay(:), '.b')
    scatter(2*ones(length(recovNight(:)),1), recovNight(:), '.k')
    errorbar([1], [mean(recovDay(:))], [std(recovDay(:))], 'ob', 'MarkerSize', 8)
    errorbar([2], [mean(recovNight(:))], [std(recovNight(:))], 'ok', 'MarkerSize', 8)
    plot([1, 2], [mean(recovDay(:)), mean(recovNight(:))], '--k', 'LineWidth', 1)
    plot([0.5, 2.5], [0,0], ':k')
    xlim([0.5, 2.5])
    ylim([-0.3,0.3]);
    ylabel('Change in Motif Frac', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', [1, 2], 'XTickLabels', [{'Over Day'}, {'Over Night'}], 'YTick', [-.3, 0, .3])
    
    %Scatterplot the recovery of Motif Continuity
    subSet = MCsnipsM((typeMask),:);
    x = xs(typeMask);
    
    h(3) = figure(3); clf
    subplot(2,2,1:2); cla
    hold on
    %Real data
    for j = 1:size(subSet,1)
        scatter(x(j).D, subSet(j, x(j).D)/subSet(j, 1), 'sb')
        scatter(x(j).N, subSet(j, x(j).N)/subSet(j, 1), 'sk')
        scatter(x(j).L, subSet(j, x(j).L)/subSet(j, 1), 'sr')      
        t(j,:) = subSet(j,:)./subSet(j, 1);
    end
    errorbar(x(j).D, mean(t(:,x(j).D)), std(t(:,x(j).D)), 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(x(j).N, mean(t(:,x(j).N)), std(t(:,x(j).N)), 'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    errorbar(x(j).L, mean(t(:,x(j).L)), std(t(:,x(j).L)), 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    plot(x(j).all, mean(t(:, x(j).all)), '--k')
    
    %Fake control data
    fakex = [1, 8, 9;1, 8, 9;1, 8, 9];
    fakey = [1, 0.92, 0.93; 1, 1.02, 1; 1, 0.93, 0.91];
    for j = 1:size(subSet,1)
        scatter(fakex(j, 1:2), fakey(j, 1:2), '.b')
        scatter(fakex(j, 3), fakey(j, 3), '.k')      
    end
    errorbar(fakex(j, 1:2), mean(fakey(:, 1:2),1), std(fakey(:, 1:2),1), 'db', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
    errorbar(fakex(j, 3), mean(fakey(:,3), 1), std(fakey(:,3), 1), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
    plot(fakex(j, :), mean(fakey, 1), ':k')
    xlim([0.5, 9.5])
    ylim([0.4,1.1]);
    ylabel('Recovery of Motif Cont', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', x(j).all, 'XTickLabels', [{'Pre'}, {'PL'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}, {'P4D'}, {'P4N'}])
    
    subplot(2,2,3); cla
    hold on
    day = 2:2:8; night = 3:2:7;
    recovDay = subSet(:, day+1) - subSet(:, day);
    recovNight = subSet(:, night+1) - subSet(:, night);
    scatter(ones(length(recovDay(:)),1), recovDay(:), '.b')
    scatter(2*ones(length(recovNight(:)),1), recovNight(:), '.k')
    errorbar([1], [mean(recovDay(:))], [std(recovDay(:))], 'ob', 'MarkerSize', 8)
    errorbar([2], [mean(recovNight(:))], [std(recovNight(:))], 'ok', 'MarkerSize', 8)
    plot([1, 2], [mean(recovDay(:)), mean(recovNight(:))], '--k', 'LineWidth', 1)
    plot([0.5, 2.5], [0,0], ':k')
    xlim([0.5, 2.5])
    ylim([-0.3,0.3]);
    ylabel('Change in Motif Cont', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', [1, 2], 'XTickLabels', [{'Over Day'}, {'Over Night'}], 'YTick', [-.3, 0, .3])

end

%Save figure and data
% saveFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
% saveName = 'Summary Electro Lesion Song Data';
% 
% %Save summary data
% save([saveFolder, saveName '.mat']);
% 
% %Save figures
% savefig(h, [saveFolder, saveName '.fig']);








