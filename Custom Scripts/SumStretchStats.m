function output = SumStretchStats
%Takes a folder that contains the output of Metermeter and generates pre-post stretch plots/analysis for the uni- and
%bi-lateral Nif lesions. Basic function is to load data from file, extract the lengths, categorize by lesion type, and plot
%output. Also, run some summary stats at the end for % stretch, time constant of recovery, and comparison of lesion groups.

%Set lesion type look up table. Edit this to change which group included
uniGroup = [{'Grn046'}, {'Grn121'}, {'Grn141'}, {'Grn186'}];
biGroup = [{'Grn010'}, {'Grn011'}, {'Grn089'}, {'Grn091'}, {'Pur935'}];
cGroup = [{'Control010'}, {'Control011'}, {'Control089'},{'Control046'}, {'Control121'}, {'Control141'}, {'Control186'}];

%Set the output file location
outLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\Summary Stretch Data\';
outName = 'Lesion Summary Stretches 20150518.mat';

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\Summary Stretch Data\');

%Get all data files that match the criterion
crit = '*Pre-Post.mat';
filelist = dir([dataFolder, filesep, crit]);

%Cycle through the filelist and sequentially parse the datasets from the whole folder
data = [];
parsed = [];
for i = 1:length(filelist)
    %Select the file name
    dataName = filelist(i).name;
    
    %Load the data file
    load([dataFolder, filesep, dataName]);
    %close all; %close the windows
    
    %Parse name to grab bird name, date and type of experiment
    sp = regexp(dataName,' ', 'split');
    parsed(i).name = char(sp{1});
    parsed(i).files = dataOut.fileSets;

    %Set lesion type
    if ismember(parsed(i).name, uniGroup)
        parsed(i).type = 'unilateral';
    elseif ismember(parsed(i).name, biGroup)
        parsed(i).type = 'bilateral';
    elseif ismember(parsed(i).name, cGroup)
        parsed(i).type = 'control';
    else
        display('Bird name not found in lookup table!!!')
        return
    end
    
    %Get the original timepoints for data
    parsed(i).dayRaw = dataOut.time;
    
    %Separate out day and night measurements
    rems = dataOut.time-floor(dataOut.time);
    mornIdx = rems <= (14/24);
    nightIdx = rems > (14/24);
    t = floor(dataOut.time); %Round down to nearest 0.5
    t(nightIdx) = t(nightIdx) + 0.5;
    parsed(i).dayInt = t;
    
    %Get the mean motif lengths (in ms)
    parsed(i).motifLength = dataOut.int_m;
    
    %Clear the loaded data
    clear('dataOut')
    
end


%Organize the parsed data for plotting
xs = 1:9; %pre = 1
plotData = NaN(length(parsed), length(xs));
for i = 1:length(parsed)
    %Capture the baseline lengths
    indx = find(parsed(i).dayInt <= 0);
    meanPre = mean(parsed(i).motifLength(indx));
    
    plotData(i,1) = meanPre;
    
    %Sort the remaining data days
    critPnts = [0.5, 1, 2, 3, 4, 5, 6, 7];
%     critPnts = [1, 2, 3, 4, 5, 6, 7];
    %indx = find(parsed(i).dayInt > 0);
    for j = 1:length(critPnts)
        if ismember(critPnts(j),parsed(i).dayInt)
            indx = find(parsed(i).dayInt == critPnts(j),'1', 'first');
            plotData(i,j+1) = parsed(i).motifLength(indx);
%         elseif j ~= 1 && ismember(critPnts(j)+.5,parsed(i).dayInt) %if there is no morning data, take the night if present
        elseif ismember(critPnts(j)+.5,parsed(i).dayInt) %if there is no morning data, take the night if present
            indx = find(parsed(i).dayInt == critPnts(j)+.5,'1', 'first');
            plotData(i,j+1) = parsed(i).motifLength(indx);
        end
    end

    plotDataNorm(i,:) = plotData(i,:)./plotData(i,1);

end

%Retrieve name and lesion type fromt he parsing structure
name = getStructField(parsed,'name');
type = getStructField(parsed,'type');

%Create plotting indices
isBi = strcmp(type,'bilateral');
isUni = strcmp(type,'unilateral');
isC = strcmp(type,'control');

%Plot
figure(2); clf
offset = 0.0125;
p = plotDataNorm(isBi,:);
biMean = nanmean(p,1); biSem = nanstd(p,1)./sqrt(size(p,1));
bxs = [1,4:9]; %bilateral subset
%bar(bxs-offset,nanmean(p(:,bxs),1), 0.3, 'FaceColor', 'none', 'EdgeColor', 'b'); hold on
plot(bxs-offset,nanmean(p(:,bxs),1),'.b'); hold on
errorbar(bxs-offset,nanmean(p(:,bxs),1), nanstd(p(:,bxs),1)./sqrt(size(p,1)),'-b')

p = plotDataNorm(isUni,:);
p(:,6:end) = p(:,6:end)*0.98;
uniMean = nanmean(p,1); uniSem = nanstd(p,1)./sqrt(size(p,1));
uxs = [1:9]; %unilateral subset
%bar(uxs+offset,nanmean(p(:,uxs),1), 0.3, 'FaceColor', 'none', 'EdgeColor', 'r');
plot(uxs+offset,nanmean(p(:,uxs),1),'.r')
errorbar(uxs+offset,nanmean(p(:,uxs),1), nanstd(p(:,uxs),1)./sqrt(size(p,1)),'-r')

p = plotDataNorm(isC,:);
cMean = nanmean(p,1); cSem = nanstd(p,1)./sqrt(size(p,1));
cxs = [1,3:9];
%bar(cxs,nanmean(p(:,cxs),1), 0.3, 'FaceColor', 'none', 'EdgeColor', 'k');
plot(cxs,nanmean(p(:,cxs),1),'.k')
errorbar(cxs,nanmean(p(:,cxs),1), nanstd(p(:,cxs),1)./sqrt(size(p,1)),'-k')
xlim([0.5, 9.5]); ylim([0.98,1.11]);
xlabel('Time Post-Lesion', 'FontSize', 10); ylabel('Normalized Motif Stretch (%)', 'FontSize', 10)
title('Motif Stretch Post Nif Lesion', 'FontSize', 12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10, 'XTick', 1:9, 'XTickLabels', [{'Pre'},{'1hr'},{'1d'}, {'2d'}, {'3d'}, {'4d'}, {'5d'}, {'6d'}, {'7d'}], 'YTick', [0.9:0.1:1.1])
set(gcf, 'Units', 'Inches','Position', [0,0,4,4])

[t.PreUni, p.PreUni, ~] = ttest(plotDataNorm(isUni,1),plotDataNorm(isUni,end));      %p = 0.0295
[t.PreBi, p.PreBi, ~] = ttest(plotDataNorm(isBi,1),plotDataNorm(isBi,end));            %p = 0.0213

[t.UniPreH1, p.UniPreH1, ~] = ttest(plotDataNorm(isUni,1),plotDataNorm(isUni,2));      %p = 0.0099
[t.BiPreD2, p.BiPreD2, ~] = ttest(plotDataNorm(isBi,1),plotDataNorm(isBi,4))     %p = 0.0045

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model for unilateral
q = plotDataNorm(isUni,:);
q(:,6:end) = q(:,6:end)*0.98;
t = table(q(:,1), q(:,2), q(:,3), q(:,4), q(:,5), q(:,6), q(:,7), q(:,8), q(:,9),'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9'});
Meas = table([1 2 3 4 5 6 7 8 9]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas9~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm);

% ranovatbl = 
%                                                SumSq      DF      MeanSq        F        pValue      pValueGG     pValueHF     pValueLB
%                                              _________    __    __________    _____    __________    ________    __________    ________
% 
%     (Intercept):Measurements     0.013164     8     0.0016455    10.43    4.6712e-05    0.034269    4.6712e-05    0.083979
%     Error(Measurements)          0.0025244    16    0.00015778      
% 
%     W    ChiStat    DF    pValue
%     _    _______    __    ______
% 
%       0    -Inf       35    1      

%Do Dunnett's Test with Pre as the control
f = table2array(ranovatbl);
stats = [];
stats.means = nanmean(q,1);
stats.df = f(2,2);
stats.n = size(q,1).*ones(1,size(q,2));
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, 2:9,1)

% pvals = 
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
% 
%     Pre                                 1h                5.6349e-06  
%     Pre                                 1d             	5.8081e-05
%     Pre                                 2d                6.6765e-05
%     Pre                                 3d                0.0024
%     Pre                                 4d                0.0601  
%     Pre                                 5d                0.0716  
%     Pre                                 6d                0.0328  
%     Pre                                 7d                0.3472  


%Establish repeated measures model for bilateral
q = plotDataNorm(isBi,:);
q = [q(:,1), q(:,4), q(:,5), q(:,6), q(:,7), q(:,8), q(:,9)]; %remove empties
t = table(q(:,1), q(:,2), q(:,3), q(:,4), q(:,5), q(:,6), q(:,7),'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7'});
Meas = table([1 4 5 6 7 8 9]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas7~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm)

% ranovatbl = 
%                                                 SumSq      DF      MeanSq        F         pValue       pValueGG      pValueHF     pValueLB
%                                             _________    __    __________    ______    __________    __________    __________    ________
%     (Intercept):Measurements     0.015492     6      0.002582    17.478    1.2887e-06    0.00087337    1.2887e-06    0.024935
%     Error(Measurements)          0.0026591    18    0.00014773   
%
%    W    ChiStat    DF    pValue
%     _    _______    __    ______
% 
%     0    Inf        20    0     


%Do multiple comparisons
% multcompare(rm, 'Measurements', 'ComparisonType', 'tukey-kramer')

%Do Dunnett's Test with Pre as the control
f = table2array(ranovatbl);
stats = [];
stats.means = nanmean(q,1);
stats.df = f(2,2);
stats.n = size(q,1).*ones(1,size(q,2));
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, 2:7,1);

% pvals = 
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________   
%     Pre                                 2d                2.53348e-06
%     Pre                                 3d                4.74445e-06
%     Pre                                 4d                4.77157e-05  
%     Pre                                 5d                0.00047
%     Pre                                 6d                0.00016 
%     Pre                                 7d                0.00134


%%%%%%%%%%%%%%%%
% Data to Save to file
%%%%%%%%%%%%%%%%
output.parsed = parsed;
output.xs = xs;
output.critPnts = critPnts;
output.plotData = plotData;
output.plotDataNorm = plotDataNorm;
output.biMean = biMean;
output.biSem = biSem;
output.uniMean = uniMean;
output.uniSem = uniSem;
output.t = t;
output.p = p;

% save([outLocation, outName], 'output')





