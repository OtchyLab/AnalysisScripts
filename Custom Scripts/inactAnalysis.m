%9/6/2014
%This script is for analying changes in song before, during, and after inactivations

%%
clear all

%Load data sources
annotName = 'V:\Grn115\Grn115_140801_annotation.mat';
fileSource = 'V:\Grn115\2014-08-01\';
inactStartTime = '12:01';
outputDir = 'C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn115\';
sp = regexp(annotName, '\', 'split');
sp2 = regexp(sp{end}, '_', 'split');
outputName = [char(sp2(1)) '_' char(sp2(2))];

load(annotName)

%%
%Define main vars
fs = 44150;
sylAbsStartTime = [];
sylDur = [];
sylType = [];
wavName = [];
sylFeat = [];
numRends = length(keys);
for i = 1:numRends
    %Unpack the time, sylType, and duration
    sylAbsStartTime = [sylAbsStartTime; elements{i}.segAbsStartTimes'];
    sylDur = [sylDur; [elements{i}.segFileEndTimes - elements{i}.segFileStartTimes]'];
    sylType = [sylType; elements{i}.segType];

    %Capture the filename of all syllables
    numSyls = length(elements{i}.segType);
    ins = cell(1,numSyls);
    for j = 1:numSyls
        ins{j} = keys{i};
    end
    wavName = [wavName; ins'];
    
    %Unpack the acoustic features...
    bDoAcoustic = true;
    if (bDoAcoustic)
        %Calc the SAP features
         features = koenigSpectral(audioread([fileSource keys{i}]), fs);
        
%         [chans, ~] = getChannels([fileSource keys{i}]);
%         features = koenigSpectral(chans(1,:), fs);

        %Parse by syllables
        ins = zeros(1,numSyls);
        for j = 1:numSyls
            s = max([1, elements{i}.segFileStartTimes(j)*1000]);
            e = min([length(features.Entropy) elements{i}.segFileEndTimes(j)*1000]);
            snip = floor(s):ceil(e);
            ins(j) = mean(features.Entropy(snip));
        end
        sylFeat = [sylFeat; ins'];
    end
    

end
sylDur = sylDur*1000; %Convert to ms
    
%%
%Define main vars
% inactStartTime = '13:07';
inactDur = 1.4; %in hours
postWash = 4; %in hours

sp = regexp(annotName, '_', 'split');
inactStartDay = sp{2};
inactDatenum = datenum([inactStartDay '_' inactStartTime], 'yymmdd_HH:MM');

%Create indices
sylInx = (sylType>=1 & sylType<101) | sylType==103;
subInx = sylType == 100;
preInx = sylAbsStartTime<inactDatenum;
% postInx = sylAbsStartTime>(inactDatenum+(inactDur/24));
postInx = sylAbsStartTime>(inactDatenum+((postWash)/24));
injInx = sylAbsStartTime>=inactDatenum & sylAbsStartTime<=(inactDatenum+(inactDur/24));

%%
%Plot duration distributions on the same axis
minval = min(sylDur);
maxval = max(sylDur);
bins = 100;
%[allX,allY] = epdf(sylDur,bins,minval,maxval);
[preX,preY] = epdf(sylDur(preInx),bins,minval,maxval);
[postX,postY] = epdf(sylDur(postInx),bins,minval,maxval);
%[subX,subY] = epdf(sylDur(subInx),bins,minval,maxval);
[injX,injY] = epdf(sylDur(injInx),bins,minval,maxval);

h(1) = figure(1);
cla; hold on
plot(preX,preY, 'LineWidth', 2)
plot(postX,postY, 'LineWidth', 2)
%plot(subX,subY, 'LineWidth', 2)
plot(injX,injY, 'LineWidth', 2)
axis tight
xlim([0 400])
ylim([0,max([preY, postY, injY])*1.1])
box off
set(gca, 'TickDir', 'out')
xlabel('Syl Dur (s)', 'FontSize', 20)
ylabel('P(t)', 'FontSize', 20)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 5])
set(gca, 'LineWidth', 3, 'FontSize', 20)
title([sp{1} ' ' sp{2}], 'FontSize', 20);
legend({'pre'; 'post'; 'inact'}); legend('boxoff')

%%
%Plot scatter for the durations
h(2) = figure(2);
cla
scatter((sylAbsStartTime-inactDatenum)*24, sylDur, '.')
hold on
line([0, 0], [0, 400], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r')
xsl = xlim;
patch([0, inactDur, inactDur, 0], [0, 0, 400, 400], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([postWash, xsl(2), xsl(2), postWash], [0, 0, 400, 400], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([xsl(1), 0, 0, xsl(1)], [0, 0, 400, 400], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
axis tight
ylim([0, 400])
box off
set(gca, 'TickDir', 'out')
xlabel('Time (hrs)', 'FontSize', 20)
ylabel('Syl Dur (ms)', 'FontSize', 20)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 7 4])
set(gca, 'LineWidth', 3, 'FontSize', 20)
title([sp{1} ' ' sp{2}], 'FontSize', 20);

%%
%Plot heatmaps for the durations and feature stats
if bDoAcoustic
    h(3) = figure(3);
    cla

    subplot(1,3,1)
    preN = ndhist(sylFeat(preInx), sylDur(preInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
    box off
    set(gca, 'TickDir', 'out')
    ylabel('Syl Dur (ms)', 'FontSize', 20)
    set(gca, 'LineWidth', 3, 'FontSize', 20)
    legend('pre'); legend('boxoff')

    subplot(1,3,2)
    injN = ndhist(sylFeat(injInx), sylDur(injInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
    box off
    set(gca, 'TickDir', 'out', 'YTickLabels', [])
    xlabel('Mean log(Entr)', 'FontSize', 20)
    set(gca, 'LineWidth', 3, 'FontSize', 20)
    title([sp{1} ' ' sp{2}], 'FontSize', 20);
    legend('elev'); legend('boxoff')

    subplot(1,3,3)
    postN = ndhist(sylFeat(postInx), sylDur(postInx), 'bins', 6, 'filter', 'axis', [-3.5 0 0 400]);
    box off
    set(gca, 'TickDir', 'out', 'YTickLabels', [])
    set(gca, 'LineWidth', 3, 'FontSize', 20)
    legend('post'); legend('boxoff')

    set(gcf, 'Units', 'Inches');
    set(gcf, 'Position', [0 0 8 5])

end
%%
%Plot heatmap of durations
% Filter parameters
window = 7;
factor = 6;

%Other parameters
maxInterval = 400;
trialBinSize = 50;
intertapBinSize = 5;
upperScale = 2.5E-2;

TargetDur = sylDur;

intertapMatrix = zeros(length(TargetDur), maxInterval);
convMatrix = ones(trialBinSize, intertapBinSize);
for m = 1:length(TargetDur)
    if (TargetDur(m) > 0) && (TargetDur(m) < maxInterval)
        intertapMatrix(m,round(TargetDur(m))) = 1;
    end
end
matrix = conv2(intertapMatrix, convMatrix, 'same');

% Filter
hfilt = ones(window*factor, window)/(factor*window^2);
matrix = imfilter(matrix, hfilt, 'replicate');
    
% Normalize
for m = 1:size(matrix,1)
    matrix(m,:) = matrix(m,:)/sum(matrix(m,:));
end

injMarkerPnt = find(sylAbsStartTime<=inactDatenum, 1, 'last');

h(4) = figure(4);
cla
imagesc(matrix, [0 upperScale]);
hold on
line([0, 400], [injMarkerPnt, injMarkerPnt], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r')

 ts = (sylAbsStartTime-inactDatenum)*24;
 Stamps = min(ts):0.25:max(ts);
for i = 1:length(Stamps)
    idxes(i) = find(ts>=Stamps(i), 1, 'first');
%     plot([(403+i), (403+i+1)],[idxes(i), idxes(i)], 'k', 'LineWidth', 4, 'clip', 'off')
%     plot((403+i),idxes(i), '.k', 'MarkerSize', 12, 'clip', 'off')
end

box off
set(gca, 'TickDir', 'out')
xlabel('Syllable Durations (ms)', 'FontSize', 20)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [0 0 5 6])
set(gca, 'LineWidth', 3, 'FontSize', 20)
title([sp{1} ' ' sp{2}], 'FontSize', 20);

%%

%Save all data
save([outputDir, outputName, '_inactDataset.mat']);

%Save figures
savefig(h, [outputDir, outputName, '.fig']);
