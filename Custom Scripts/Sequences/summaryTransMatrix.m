%Sequentially calculates transititon probability matrices (by calling
%makeTransMatrix.m). Saves summary data to file for later processing.

%% Setup the script
%clear

%Folder to load annotation files from

%ChABC in HVC
% mother = '/Users/Tim/Dropbox/Summary Ints/HVC';
% type = 'HVC_ChABC';

%PBS in HVC
% mother = '/Users/Tim/Dropbox/Summary Ints/PBS';
% type = 'HVC_PBS';

%ChABC into RA
mother = '/Users/Tim/Dropbox/Summary Ints/RA';
type = 'RA_ChABC';

annots = dir([mother filesep '*annotation.mat']);

%Where to save the output
mos = '/Users/Tim/Dropbox/Summary Ints';
outFile = 'seqStats_all_v1.mat';

%% Sequentially process the annotations

for i = 1:numel(annots)
    %Retrieve the matrices
    fname = [mother, filesep, annots(i).name];
    [seqMatrix, seqProbMatrix, sylTypes, mainTrans, figHand(end+1)] = makeTransMatrix(fname);
    
    %Parse name
    sp = regexp(annots(i).name, '_', 'split');
    
    %Populate the carried variable
    seqStats = [];
    seqStats.bird = sp{1};
    seqStats.date = sp{2};
    seqStats.type = type;
    seqStats.filename = annots(i).name;
    
    seqStats.sylTypes = sylTypes;
    seqStats.mainTrans = mainTrans;
    seqStats.seqMatrix = seqMatrix;
    seqStats.seqProbMatrix = seqProbMatrix;
    
    %Save to file
    outName = [mos, filesep, outFile];
    m = exist(outName);
    if m == 2
        %File already exists
        load(outName, 'stats')
        stats(end+1) = seqStats;
    else
        %No file yet created
        stats = seqStats;
    end
    
    %Save the updated data to file
    save(outName, 'stats')
    
    disp(['Done with ' annots(i).name])
end

%% Calculate Difference Matrix and 
diffStats = [];

bVect = getFieldVectorCell(stats, 'bird');
birdlist = unique(bVect);
bool = strcmp(birdlist, 'LW58');
birdlist(bool) = [];

for i = 1:numel(birdlist)
    %Masking
    mask = ismember(bVect, birdlist(i));
    subset = stats(mask);
    
    %makeDiff
    seqDiffMatrix = abs(subset(1).seqProbMatrix - subset(2).seqProbMatrix);
    
    %Retrieve the primary motif differences
    diffs = [];
    ts = subset(1).mainTrans;
    for j = 1:size(ts,1)
        diffs = [diffs, seqDiffMatrix(ts(j,1),ts(j,2))];
    end
    
    %Difference score as the mean across primary motif
    diffScore = mean(diffs);

    %Populate the carried variable

    diffStats(i).bird = subset(1).bird;
    diffStats(i).date = [subset(1).date, ' -> ', subset(2).date];
    diffStats(i).type = subset(1).type;
    
    diffStats(i).sylTypes = subset(1).sylTypes;
    diffStats(i).mainTrans = subset(1).mainTrans;
    diffStats(i).seqDiffMatrix = seqDiffMatrix;
    diffStats(i).diffScore = diffScore;

    %Add the difference matrix to the plot
    h = figure('Name', [subset(2).filename(1:end-15), ' DiffMat']);
    
    %Plot the counts matrix
    imagesc([1 numel(sylTypes)],[1 numel(sylTypes)], seqDiffMatrix, [0, 1]);
    colormap(bone)
    colorbar
    
    axis square
    xlabel('To')
    ylabel('From')
    set(gca, 'Box', 'off', 'TickDir', 'out')
    set(gca, 'XTick', 1:numel(sylTypes), 'XTickLabel', sylTypes)
    set(gca, 'YTick', 1:numel(sylTypes), 'YTickLabel', sylTypes)
    
    savefig(h, [mos, filesep, subset(2).filename(1:end-15), '.fig'])
end

tVect = getFieldVectorCell(diffStats, 'type');
[typelist, I] = sort(tVect);
diffStats = diffStats(I);

%Save to file
outName = [mos, filesep, 'diffStats.mat'];

%Save the updated data to file
save(outName, 'diffStats')

disp(['Done with diffStats.mat'])

%% Plot the difference scores by category

%Get indices
tVect = getFieldVectorCell(diffStats, 'type');
typelist = unique(tVect);
typelist = typelist([1, 3, 2]);
diffScoreVect = getFieldVector(diffStats, 'diffScore');
xs = 1:3;

ms = [];
ss = [];
xss = []; yss = [];
for i = 1:numel(typelist)
    %Masking
    mask = ismember(tVect, typelist(i));
    subset = diffScoreVect(mask);
    
    %Descriptive stats
    ms(i) = mean(subset);
    ss(i) = std(subset);
    
    %Scatter totals
    xss = [xss, xs(i)*ones(size(subset))];
    yss = [yss, subset];
end

%Add the difference matrix to the plot
h = figure('Name', 'Summary Difference Scores'); clf
set(gcf, 'Units', 'Inches', 'Position', [7, 6.75, 5, 3.75])

%Plot the bar
b = bar(xs, ms, 'k'); hold on
b.LineWidth = 1.5;
b.FaceColor = [0.5, 0.5, 0.5];

e = errorbar(xs, ms, ss, 'Color', 'k', 'LineStyle', 'none');
e.LineWidth = 1.5;

s = scatter(xss, yss, 'ok', 'jitter','on', 'jitterAmount',0.05);
s.SizeData = 100;

xlim([0, 4])
ylabel('Sequence Difference Score')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', xs, 'XTickLabel', typelist)
set(gca, 'YTick', [0, 0.5])

% savefig(h, [mos, filesep, subset(2).filename(1:end-15), '.fig'])


