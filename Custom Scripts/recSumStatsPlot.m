function recSumStatsPlot(longRecStats)
 close all
summary = [];

name = getStructField(longRecStats,'name');
type = getStructField(longRecStats,'type');
xs = getStructField(longRecStats,'xs');

meanPowStatM = getStructField(longRecStats,'meanPowStatM');
meanPowStatS = getStructField(longRecStats,'meanPowStatS');
WinCorrM = getStructField(longRecStats,'WinCorrM');
WinCorrS = getStructField(longRecStats,'WinCorrS');
neuroCorr = getStructField(longRecStats,'neuroCorr');

labels =  [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}];

birds = unique(name);
conditions = [{'hit'}, {'control'}];
colors = [{'b'}, {'r'}];

for i = 1:length(conditions)
    typeMask = strcmp(conditions(i), type);
    
    %Scatterplot the recovery of mean power
    subSet = meanPowStatM((typeMask),:);
    semNum = sqrt(size(subSet,1));
    x = xs(typeMask);
    
    h(1) = figure(1);
    subplot(2,2,1:2);
    hold on    
    %Real data
    P = [];
    for j = 1:size(subSet,1)
        plot(x(j).all, subSet(j, x(j).all)/subSet(j, 1), 'o', 'Color', colors{i}, 'MarkerSize', 6)
        P(j,:) = subSet(j,:)./subSet(j, 1);
    end
    errorbar(x(j).all, mean(P(:,x(j).all)), std(P(:,x(j).all))./semNum, colors{i}, 'MarkerSize', 8)
    plot(x(j).all, mean(P(:, x(j).all)), ['--', colors{i}])
    xlim([0.5, 9.5])
    ylim([0.2,1.1]);
    ylabel('Norm. HVC Activity Power', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', x(j).all, 'XTickLabels', [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}])
    
    if i ==1
        P_lesion = P;
        subplot(2,2,3); cla
        hold on
        recovPowDay = []; recovPowNight = [];
        for j = 1:size(subSet,1)
            %Mean change over 3 days or 2 nights (normalized
            %Day
            recovPowDay(j,1) = (subSet(j, 3) - subSet(j, 2))./subSet(j, 9);
            recovPowDay(j,2) = (subSet(j, 5) - subSet(j, 4))./subSet(j, 9);
            recovPowDay(j,3) = (subSet(j, 7) - subSet(j, 6))./subSet(j, 9);
            
            %Night
            recovPowNight(j,1) = (subSet(j, 4) - subSet(j, 3))./subSet(j, 9);
            recovPowNight(j,2) = (subSet(j, 6) - subSet(j, 5))./subSet(j, 9);
        end
        
        plot(ones(j,1), mean(recovPowDay,2), 'o', 'Color', [0.5, 0.5, 0.5],  'MarkerSize', 6)
        plot(2*ones(j,1), mean(recovPowNight,2), 'o', 'Color', [0.5, 0.5, 0.5],  'MarkerSize', 6)
        errorbar([1, 2], [mean(mean(recovPowDay,2)), mean(mean(recovPowNight,2))], [std(mean(recovPowDay,2)), std(mean(recovPowNight,2))]./semNum, '.k', 'MarkerSize', 10)
        bar([1, 2], [mean(mean(recovPowDay,2)), mean(mean(recovPowNight,2))], 'FaceColor', 'none')
%         plot([0.5, 2.5], [0,0], ':k')
        xlim([0.5, 2.5])
        ylim([-0.1,0.25]);
        ylabel('Change in Norm Power', 'FontSize', 10)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
        set(gca, 'XTick', [1, 2], 'XTickLabels', [{'Over Day'}, {'Over Night'}], 'YTick', [-0.1:0.1:0.3])
    end
    
    set(gcf, 'Units', 'Inches', 'Position', [0 0 4 6])
    
    %Scatterplot the recovery of Neural dynamics
    subSet = neuroCorr((typeMask),:);
    x = xs(typeMask);
    h(2) = figure(2);
    subplot(2,2,1:2);
    hold on
    C = [];
    %Real data
    for j = 1:size(subSet,1)

        %Real
        plot(x(j).all, subSet(j, x(j).all), 'o', 'Color', colors{i}, 'MarkerSize', 6)
        C(j,:) = subSet(j,:);
    end
    errorbar(x(j).all, mean(C(:,x(j).all)), std(C(:,x(j).all))./semNum, colors{i}, 'MarkerSize', 8)
    plot(x(j).all, mean(C(:, x(j).all)), ['--' colors{i}])
    
    %Format
    xlim([0.5, 9.5])
    ylim([0.1,1.1]);
    ylabel('Correlation to Pre-Lesion', 'FontSize', 10)
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
    set(gca, 'XTick', x(j).all, 'XTickLabels', [{'Pre'}, {'PL'}, {'PN'}, {'P1D'}, {'P1N'}, {'P2D'}, {'P2N'}, {'P3D'}, {'P3N'}])
    
    if i ==1
        C_lesion = C;
        subplot(2,2,3); cla
        hold on
        recovCorrDay = []; recovCorrNight = [];
        for j = 1:size(subSet,1)
            %Mean change over 3 days or 2 nights (normalized
            %Day
            recovCorrDay(j,1) = (subSet(j, 3) - subSet(j, 2))./subSet(j, 9);
            recovCorrDay(j,2) = (subSet(j, 5) - subSet(j, 4))./subSet(j, 9);
            recovCorrDay(j,3) = (subSet(j, 7) - subSet(j, 6))./subSet(j, 9);
            
            %Night
            recovCorrNight(j,1) = (subSet(j, 4) - subSet(j, 3))./subSet(j, 9);
            recovCorrNight(j,2) = (subSet(j, 6) - subSet(j, 5))./subSet(j, 9);
        end
        plot(ones(j,1), mean(recovCorrDay,2), 'o', 'Color', [0.5, 0.5, 0.5],  'MarkerSize', 6)
        plot(2*ones(j,1), mean(recovCorrNight,2), 'o', 'Color', [0.5, 0.5, 0.5],  'MarkerSize', 6)
        errorbar([1, 2], [mean(mean(recovCorrDay,2)), mean(mean(recovCorrNight,2))], [std(mean(recovCorrDay,2)), std(mean(recovCorrNight,2))]./semNum, '.k', 'MarkerSize', 10)
        bar([1, 2], [mean(mean(recovCorrDay,2)), mean(mean(recovCorrNight,2))], 'FaceColor', 'none')
        
        %Format
        xlim([0.5, 2.5])
        ylim([0, 0.3])
        ylabel('Change in X-Rend Corr', 'FontSize', 10)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10)
        set(gca, 'XTick', [1, 2], 'XTickLabels', [{'Over Day'}, {'Over Night'}], 'YTick', [0:0.1:0.3])
    end

        set(gcf, 'Units', 'Inches', 'Position', [0 0 4 6])
        
end

%Hypothesis testing
[t.corrDayNight,p.corrDayNight,~] = ttest(mean(recovCorrDay,2),mean(recovCorrNight,2));         %p = 0.0226
[t.powDayNight,p.powDayNight,~] = ttest(mean(recovPowDay,2),mean(recovPowNight,2));        %p = 0.0302
[t.corrDay,p.corrDay,~] = ttest(mean(recovCorrDay,2));                                                          %p = 0.0532
[t.corrNight,p.corrNight,~] = ttest(mean(recovCorrNight,2));                                                     %p = 0.0098
[t.powDay,p.powDay,~] = ttest(mean(recovPowDay,2));                                                         %p = 0.0070
[t.powNight,p.powNight,~] = ttest(mean(recovPowNight,2));                                                    %p = 0.4787

%Corr recovery
% [t.corrPost1h,p.corrPost1h,~] = ttest2(P_lesion(:,2),P(:,2));     %p = 0.0039
% [t.corrPost8h,p.corrPost8h,~] = ttest2(P_lesion(:,3),P(:,3));     %p = 0.0017
% [t.corrPost1d,p.corrPost1d,~] = ttest2(P_lesion(:,4),P(:,4));     %p = 0.0012
% [t.corrPost1n,p.corrPost1n,~] = ttest2(P_lesion(:,5),P(:,5));     %p = 6.493e-4
% [t.corrPost2d,p.corrPost2d,~] = ttest2(P_lesion(:,6),P(:,6));     %p = 0.0300
% [t.corrPost2n,p.corrPost2n,~] = ttest2(P_lesion(:,7),P(:,7));     %p = 0.0107
% [t.corrPost3d,p.corrPost3d,~] = ttest2(P_lesion(:,8),P(:,8));     %p = 0.0969
% [t.corrPost3n,p.corrPost3n,~] = ttest2(P_lesion(:,9),P(:,9));     %p = 0.6353
% 
% [t.powPost1h,p.powPost1h,~] = ttest2(C_lesion(:,2),C(:,2));     %p = 6.0037e-4
% [t.powPost8h,p.powPost8h,~] = ttest2(C_lesion(:,3),C(:,3));     %p = 0.0019
% [t.powPost1d,p.powPost1d,~] = ttest2(C_lesion(:,4),C(:,4));     %p = 0.0180
% [t.powPost1n,p.powPost1n,~] = ttest2(C_lesion(:,5),C(:,5));     %p = 0.0309
% [t.powPost2d,p.powPost2d,~] = ttest2(C_lesion(:,6),C(:,6));     %p = 0.0165
% [t.powPost2n,p.powPost2n,~] = ttest2(C_lesion(:,7),C(:,7));     %p = 0.0073
% [t.powPost3d,p.powPost3d,~] = ttest2(C_lesion(:,8),C(:,8));     %p = 0.0255
% [t.powPost3n,p.powPost3n,~] = ttest2(C_lesion(:,9),C(:,9))      %p = 0.2857

[t.corrPost1h,p.corrPost1h,~] = ttest(P_lesion(:,2),P_lesion(:,1));     %p = 0.0031
[t.corrPost8h,p.corrPost8h,~] = ttest(P_lesion(:,3),P_lesion(:,1));     %p = 0.0031
[t.corrPost1d,p.corrPost1d,~] = ttest(P_lesion(:,4),P_lesion(:,1));     %p = 0.0028
[t.corrPost1n,p.corrPost1n,~] = ttest(P_lesion(:,5),P_lesion(:,1));     %p = 0.0028
[t.corrPost2d,p.corrPost2d,~] = ttest(P_lesion(:,6),P_lesion(:,1));     %p = 0.0122
[t.corrPost2n,p.corrPost2n,~] = ttest(P_lesion(:,7),P_lesion(:,1));     %p = 0.0142
[t.corrPost3d,p.corrPost3d,~] = ttest(P_lesion(:,8),P_lesion(:,1));     %p = 0.0099
[t.corrPost3n,p.corrPost3n,~] = ttest(P_lesion(:,9),P_lesion(:,1));     %p = 0.3339

[t.powPost1h,p.powPost1h,~] = ttest(C_lesion(:,2),C_lesion(:,1));     %p = 0.0022
[t.powPost8h,p.powPost8h,~] = ttest(C_lesion(:,3),C_lesion(:,1));     %p = 0.0043
[t.powPost1d,p.powPost1d,~] = ttest(C_lesion(:,4),C_lesion(:,1));     %p = 0.0173
[t.powPost1n,p.powPost1n,~] = ttest(C_lesion(:,5),C_lesion(:,1));     %p = 0.0157
[t.powPost2d,p.powPost2d,~] = ttest(C_lesion(:,6),C_lesion(:,1));     %p = 0.0067
[t.powPost2n,p.powPost2n,~] = ttest(C_lesion(:,7),C_lesion(:,1));     %p = 0.0026
[t.powPost3d,p.powPost3d,~] = ttest(C_lesion(:,8),C_lesion(:,1));     %p = 0.0084
[t.powPost3n,p.powPost3n,~] = ttest(C_lesion(:,9),C_lesion(:,1));      %p = 0.0099


% %Run one way, repeated measures ANOVA
% %Establish repeated measures model
% t = table(ddDist, inactDist, PBSDist, elevDist,ppDist,'VariableNames',{'cond1', 'cond2', 'cond3', 'cond4', 'cond5'});
% Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
% rm = fitrm(t,'cond1-cond5~1','WithinDesign',Meas);
% 
% %Do the RM ANOVA
% ranovatbl = ranova(rm)
% f = table2array(ranovatbl);
% stats = [];
% stats.means = mean([ddDist, inactDist, PBSDist, elevDist,ppDist],1);
% stats.df = f(2,2);
% stats.n = [size(ddDist,1), size(inactDist,1), size(PBSDist,1), size(elevDist,1), size(ppDist,1),];
% stats.s = sqrt(f(2,3));
































%Save figure and data
% saveFolder = 'C:\Users\Tim\Desktop\Nif Project Figures\ElectroLesions\';
% saveName = 'Summary Electro Lesion Rec Data V2';

% %Save summary data
% save([saveFolder, saveName '.mat']);
% 
% %Save figures
% savefig(h, [saveFolder, saveName '.fig']);







