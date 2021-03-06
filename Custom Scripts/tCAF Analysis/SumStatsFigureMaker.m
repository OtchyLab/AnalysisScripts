%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes timeseries plots of the syllable durations
%and changes. Also saves some summary data to a file to group data over
%birds and conditions.

%% Set it up
%clear

%File location
mother = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis';
% mother = '/Users/Tim/Documents/MATLAB/General/Custom Scripts/tCAF Analysis';
file = 'cleanedSummaryStatsv2.mat'; %<=== Update this to load different file

%Load it
%load([mother, filesep, file])

%Plot List
cats = {'Control', 'ChABC'};
mkrs = {'o', 's'};
stretchCents =[1, 2, 4, 5];
rateCents =[1, 2, 4, 5, 7, 8];
col = {'b', 'r'};

%% Create the figure
f = figure(102); clf
set(f, 'Units', 'inches', 'Position', [0.5, 3.5, 6.5, 7])

%Plot max streches
subplot(2,1,1); cla

%Group data
numCats = numel(cats);
v = getFieldVectorCell(stats, 'batch');
for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    
    %apply
    subset = stats(mask);
    numSets = numel(subset);
    tStretch = []; nStretch = [];
    for j = 1:numSets
        %Running list
        tStretch = [tStretch,  subset(j).stretch.target'];
        nStretch = [nStretch, subset(j).stretch.nontarget(:)'];
        
        %Plot the stretches
        tXs = stretchCents(i)*ones(size(subset(j).stretch.target));
        nXs = stretchCents(i+2)*ones(size(subset(j).stretch.nontarget(:)));
        scatter([tXs', nXs'], [subset(j).stretch.target', subset(j).stretch.nontarget(:)'],'k', 'Marker', mkrs{j}); hold on
    end
    
    bar(stretchCents([i, i+2]), [mean(tStretch), mean(nStretch)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
    errorbar(stretchCents([i, i+2]), [mean(tStretch), mean(nStretch)], [std(tStretch, 1), std(nStretch,1)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
    
    
end
%Format
xlim([0, 6])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'Box', 'off', 'XTick', stretchCents, 'XTickLabel', {'Intact', 'ChABC', 'Intact', 'ChABC'})
ylabel('Peak Stretch (ms)')
xlabel('Target                                                Non-Target')

%Plot shift rates
subplot(2,1,2); cla

for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    
    %apply
    subset = stats(mask);
    numSets = numel(subset);
    dUp = []; dDown = []; dSpont = [];
    for j = 1:numSets
        %Running list
        dUp = [dUp,  subset(j).rates.up];
        dDown = [dDown, subset(j).rates.down];
        dSpont = [dSpont, subset(j).rates.spon];
        
        %Plot the stretches
        uXs = rateCents(i)*ones(size(subset(j).rates.up));
        dXs = rateCents(i+2)*ones(size(subset(j).rates.down));
        sXs = rateCents(i+4)*ones(size(subset(j).rates.spon));
        scatter([uXs, dXs, sXs], [subset(j).rates.up, subset(j).rates.down, subset(j).rates.spon],'k', 'Marker', mkrs{j}); hold on
    end
    
    bar(rateCents([i, i+2, i+4]), [mean(dUp), mean(dDown), mean(dSpont)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
    errorbar(rateCents([i, i+2, i+4]), [mean(dUp), mean(dDown), mean(dSpont)], [std(dUp), std(dDown), std(dSpont)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
    
end
%Format
xlim([0, 9])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'Box', 'off', 'XTick', rateCents, 'XTickLabel', {'Intact', 'ChABC', 'Intact', 'ChABC', 'Intact', 'ChABC'})
ylabel('Stretch rate (ms/day)')
xlabel('Drive Up                        Drive Down                        Spont Down')


%% Create the % figure
f = figure(103); clf
set(f, 'Units', 'inches', 'Position', [0.5, 3.5, 6.5, 7])

%Plot max streches
subplot(2,1,1); cla

%Group data
numCats = numel(cats);
v = getFieldVectorCell(stats, 'batch');
for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    
    %apply
    subset = stats(mask);
    numSets = numel(subset);
    tStretch = []; nStretch = [];
    for j = 1:numSets
        %Running list
        tStretch = [tStretch,  100.*(subset(j).stretch.target./subset(j).bdur.target)'];
        nStretch = [nStretch, 100.*(subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'];
        
        %Plot the stretches
        tXs = stretchCents(i)*ones(size(subset(j).stretch.target));
        nXs = stretchCents(i+2)*ones(size(subset(j).stretch.nontarget(:)));
        scatter([tXs', nXs'], 100.*[(subset(j).stretch.target./subset(j).bdur.target)', (subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'],'k', 'Marker', mkrs{j}); hold on
    end
    
    bar(stretchCents([i, i+2]), [mean(tStretch), mean(nStretch)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
    errorbar(stretchCents([i, i+2]), [mean(tStretch), mean(nStretch)], [std(tStretch, 1), std(nStretch,1)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
    
    
end
%Format
xlim([0, 6])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'Box', 'off', 'XTick', stretchCents, 'XTickLabel', {'Intact', 'ChABC', 'Intact', 'ChABC'})
ylabel('Max Stretch (%)')
xlabel('Target                                                Non-Target')

%Plot shift rates
subplot(2,1,2); cla

for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    
    %apply
    subset = stats(mask);
    numSets = numel(subset);
    dUp = []; dDown = []; dSpont = [];
    for j = 1:numSets
        %Running list
        dUp = [dUp,  100.*(subset(j).rates.up'./subset(j).bdur.target)'];
        dDown = [dDown, 100.*subset(j).rates.down./subset(j).bdur.target(1)];
        dSpont = [dSpont, 100.*subset(j).rates.spon./subset(j).bdur.target(2)];
        
        %Plot the stretches
        uXs = rateCents(i)*ones(size(subset(j).rates.up));
        dXs = rateCents(i+2)*ones(size(subset(j).rates.down));
        sXs = rateCents(i+4)*ones(size(subset(j).rates.spon));
        scatter([uXs, dXs, sXs], 100.*[(subset(j).rates.up'./subset(j).bdur.target)', subset(j).rates.down./subset(j).bdur.target(1), subset(j).rates.spon./subset(j).bdur.target(2)],'k', 'Marker', mkrs{j}); hold on
    end
    
    bar(rateCents([i, i+2, i+4]), [mean(dUp), mean(dDown), mean(dSpont)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
    errorbar(rateCents([i, i+2, i+4]), [mean(dUp), mean(dDown), mean(dSpont)], [std(dUp), std(dDown), std(dSpont)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
    
end
%Format
xlim([0, 9])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'Box', 'off', 'XTick', rateCents, 'XTickLabel', {'Intact', 'ChABC', 'Intact', 'ChABC', 'Intact', 'ChABC'})
ylabel('Stretch rate (%/day)')
xlabel('Drive Up                        Drive Down                        Spont Down')

%% Create the % breakout figure
f = figure(104); clf
set(f, 'Units', 'inches', 'Position', [6, 5.5, 5.5, 7.5])
% stretchCents =[1, 2, 4, 5];
stretchCents = [1, 2, 3, 5, 6, 7];
rateCents =[1, 2, 3, 5, 6, 7, 9, 10, 11];

%Plot max streches
subplot(2,1,1); cla

%Group data
numCats = numel(cats);
v = getFieldVectorCell(stats, 'batch');
for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    subset = stats(mask);
    numSets = numel(subset);
    
    if i == 1 %Control plots
        tStretch = []; nStretch = [];
        for j = 1:numSets
            %Running list
            tStretch = [tStretch,  100.*(subset(j).stretch.target./subset(j).bdur.target)'];
            nStretch = [nStretch, 100.*(subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'];
            
            %Plot the stretches
            tXs = stretchCents(i)*ones(size(subset(j).stretch.target));
            nXs = stretchCents(i+3)*ones(size(subset(j).stretch.nontarget(:)));
            scatter([tXs', nXs'], 100.*[(subset(j).stretch.target./subset(j).bdur.target)', (subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'],'k', 'Marker', mkrs{j}, 'jitter','on', 'jitterAmount',0.05); hold on
        end
        
%         bar(stretchCents([i, i+3]), [mean(tStretch), mean(nStretch)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
%         errorbar(stretchCents([i, i+3]), [mean(tStretch), mean(nStretch)], [std(tStretch, 1), std(nStretch,1)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
%         
    else %ChABC plots
        
        toStretch = []; tnStretch = []; naStretch = []; nbStretch = [];
        for j = 1:numSets
            
            %Running list
            if strcmp(subset(j).bird, 'LW60') %no switch
                toStretch = [toStretch,  100.*(subset(j).stretch.target./subset(j).bdur.target)'];
                naStretch = [naStretch, 100.*(subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'];
                tXs = stretchCents(i)*ones(size(subset(j).stretch.target));
                nXs = stretchCents(i+3)*ones(size(subset(j).stretch.nontarget(:)));
            elseif strcmp(subset(j).bird, 'LW58') %switch
                tnStretch = [tnStretch,  100.*(subset(j).stretch.target./subset(j).bdur.target)'];
                nbStretch = [nbStretch, 100.*(subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'];
                tXs = stretchCents(i+1)*ones(size(subset(j).stretch.target));
                nXs = stretchCents(i+4)*ones(size(subset(j).stretch.nontarget(:)));
            end

            scatter([tXs', nXs'], 100.*[(subset(j).stretch.target./subset(j).bdur.target)', (subset(j).stretch.nontarget(:)./subset(j).bdur.nontarget(:))'],'k', 'Marker', mkrs{j}, 'jitter','on', 'jitterAmount',0.05); hold on
        end
        
        bar(stretchCents, [mean(tStretch), mean(toStretch), mean(tnStretch), mean(nStretch), mean(naStretch),mean(nbStretch)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.6)
        errorbar(stretchCents, [mean(tStretch), mean(toStretch), mean(tnStretch), mean(nStretch),mean(naStretch),mean(nbStretch)], [std(tStretch, 1), std(toStretch, 1), std(tnStretch, 1),std(nStretch,1),std(naStretch,1),std(nbStretch,1)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
        
%         bar(stretchCents([i, i+1, i+3, i+4]), [mean(toStretch), mean(tnStretch),mean(naStretch),mean(nbStretch)], col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.3)
%         errorbar(stretchCents([i, i+1, i+3, i+4]), [mean(toStretch), mean(tnStretch),mean(naStretch),mean(nbStretch)], [std(toStretch, 1), std(tnStretch, 1),std(naStretch,1),std(nbStretch,1)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
%         
        
    end
    
end
%Format
xlim([0, 8])
ylim([-5, 20])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'Box', 'off', 'YTick', -5:10:25, 'XTick', stretchCents, 'XTickLabel', {'Baseline', 'Stable', 'Re-Seq'})
ylabel('Max Stretch (%)')
xlabel('Target                          Non-Target')

%Plot shift rates
subplot(2,1,2); cla

for i = 1:numCats
    %Make analysis mask
    mask = strcmp(v, cats{i});
    subset = stats(mask);
    numSets = numel(subset);
    
    %apply
    if i == 1 %Control plots
        dUp = []; dDown = []; dSpont = [];
        for j = 1:numSets
            %Running list
            dUp = [dUp,  100.*(subset(j).rates.up'./subset(j).bdur.target)'];
            dDown = [dDown, 100.*subset(j).rates.down./subset(j).bdur.target(1)];
            dSpont = [dSpont, 100.*subset(j).rates.spon./subset(j).bdur.target(2)];
            
            %Plot the stretches
            uXs = rateCents(i)*ones(size(subset(j).rates.up));
            dXs = rateCents(i+3)*ones(size(subset(j).rates.down));
            sXs = rateCents(i+6)*ones(size(subset(j).rates.spon));
            scatter([uXs, dXs, sXs], abs(100.*[(subset(j).rates.up'./subset(j).bdur.target)', subset(j).rates.down./subset(j).bdur.target(1), subset(j).rates.spon./subset(j).bdur.target(2)]),'k', 'Marker', mkrs{j}, 'jitter','on', 'jitterAmount',0.05); hold on
        end
        
        bar(rateCents([i, i+3, i+6]), abs([mean(dUp), mean(dDown), mean(dSpont)]), col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.2)
        errorbar(rateCents([i, i+3, i+6]), abs([mean(dUp), mean(dDown), mean(dSpont)]), [std(dUp), std(dDown), std(dSpont)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
        
    elseif i == 2 %ChABC plots
        for j = 1:numSets
            dUp = []; dDown = []; dSpont = [];
            %Running list
            if strcmp(subset(j).bird, 'LW60') %no switch
                dUp = [dUp,  100.*(subset(j).rates.up'./subset(j).bdur.target)'];
                dDown = [dDown, 100.*subset(j).rates.down./subset(j).bdur.target(1)];
                dSpont = [dSpont, 100.*subset(j).rates.spon./subset(j).bdur.target(2)];
                
                uXs = rateCents(i)*ones(size(subset(j).rates.up));
                dXs = rateCents(i+3)*ones(size(subset(j).rates.down));
                sXs = rateCents(i+6)*ones(size(subset(j).rates.spon));
                
                scatter([uXs, dXs, sXs], abs(100.*[(subset(j).rates.up'./subset(j).bdur.target)', subset(j).rates.down./subset(j).bdur.target(1), subset(j).rates.spon./subset(j).bdur.target(2)]),'k', 'Marker', mkrs{j}, 'jitter','on', 'jitterAmount',0.05); hold on
                bar(rateCents([i, i+3, i+6]), abs([mean(dUp), mean(dDown), mean(dSpont)]), col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.2)
        errorbar(rateCents([i, i+3, i+6]), abs([mean(dUp), mean(dDown), mean(dSpont)]), [std(dUp), std(dDown), std(dSpont)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
        
            elseif strcmp(subset(j).bird, 'LW58') %switch
                dUp = [dUp,  100.*(subset(j).rates.up'./subset(j).bdur.target)'];
                dDown = [dDown, 100.*subset(j).rates.down./subset(j).bdur.target(1)];
                dSpont = [dSpont, 100.*subset(j).rates.spon./subset(j).bdur.target(2)];
                
                uXs = rateCents(i+1)*ones(size(subset(j).rates.up));
                dXs = rateCents(i+3+1)*ones(size(subset(j).rates.down));
                sXs = rateCents(i+6+1)*ones(size(subset(j).rates.spon));
                
                scatter([uXs, dXs, sXs], abs(100.*[(subset(j).rates.up'./subset(j).bdur.target)', subset(j).rates.down./subset(j).bdur.target(1), subset(j).rates.spon./subset(j).bdur.target(2)]),'k', 'Marker', mkrs{j}, 'jitter','on', 'jitterAmount',0.05); hold on
                bar(rateCents([i+1, i+3+1, i+6+1]), abs([mean(dUp), mean(dDown), mean(dSpont)]), col{i}, 'FaceAlpha', 0.6, 'BarWidth', 0.2)
        errorbar(rateCents([i+1, i+3+1, i+6+1]), abs([mean(dUp), mean(dDown), mean(dSpont)]), [std(dUp), std(dDown), std(dSpont)], col{i}, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none')
            end
            
        end
        
        
    end
end
%Format
xlim([0, 12])
ylim([0, 3.5])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'YTick', [0,3], 'XTick', rateCents, 'XTickLabel', {'Baseline', 'Stable', 'Re-Seq'})
ylabel('Stretch rate (%/day)')
xlabel('Drive Up           Drive Down            Spont Down')












