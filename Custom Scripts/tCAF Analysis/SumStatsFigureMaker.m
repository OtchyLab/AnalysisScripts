%tCAF analysis script for Nerurkar and Otchy (2018) and Darkwa, Nerurkar,
%and Otchy (2018). This makes timeseries plots of the syllable durations
%and changes. Also saves some summary data to a file to group data over
%birds and conditions.

%% Set it up
%clear 

%File location
mother = 'C:\Users\Tim\Desktop\Matlab Code\General Scripts\Custom Scripts\tCAF Analysis';
% mother = '/Users/Tim/Documents/MATLAB/General/Custom Scripts/tCAF Analysis';
file = 'summaryStats.mat'; %<=== Update this to load different file

%Load it
%load([mother, filesep, file])

cats = {'Control', 'ChABC'};

%Create the figure
f = figure(102); clf
set(f, 'Units', 'inches', 'Position', [0.5, 3.5, 6.5, 7])

%Plot List
mkrs = {'o', 'd', '+'};
stretchCents =[1, 2, 4, 5];
rateCents =[1, 2, 4, 5, 7, 8];
col = {'b', 'r'};

%% Plot max streches
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
        tStretch = [tStretch,  subset(j).tStretch'];
        nStretch = [nStretch, subset(j).ntStretch(:)'];
        
        %Plot the stretches
        tXs = stretchCents(i)*ones(size(subset(j).tStretch));
        nXs = stretchCents(i+2)*ones(size(subset(j).ntStretch(:)));
        scatter([tXs', nXs'], [subset(j).tStretch', subset(j).ntStretch(:)'],'k', 'Marker', mkrs{j}); hold on
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


%% Plot shift rates
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
        dUp = [dUp,  subset(j).shiftRate([1, 4], end-1)'];
        dDown = [dDown, subset(j).shiftRate(2, end-1)'];
        dSpont = [dSpont, subset(j).shiftRate(5, end-1)'];
        
        %Plot the stretches
        uXs = rateCents(i)*ones(size(subset(j).shiftRate([1, 4], end-1)));
        dXs = rateCents(i+2)*ones(size(subset(j).shiftRate(2, end-1)));
        sXs = rateCents(i+4)*ones(size(subset(j).shiftRate(5, end-1)));
        scatter([uXs', dXs', sXs'], [subset(j).shiftRate([1, 4], end-1)', subset(j).shiftRate(2, end-1)', subset(j).shiftRate(5, end-1)'],'k', 'Marker', mkrs{j}); hold on
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
















