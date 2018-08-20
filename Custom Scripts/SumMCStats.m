function output = SumMCStats
%Takes a folder that contains the output of Metermeter and generates pre-post stretch plots/analysis for the uni- and
%bi-lateral Nif lesions. Basic function is to load data from file, extract the lengths, categorize by lesion type, and plot
%output. Also, run some summary stats at the end for % stretch, time constant of recovery, and comparison of lesion groups.

%Set lesion type look up table. Edit this to change which group included
uniGroup = [{'Grn046'}, {'Grn121'}, {'Grn141'}, {'Grn186'}];
biGroup = [{'Grn010'}, {'Grn011'}, {'Grn089'}, {'Grn091'}, {'Pur935'}];
cGroup = [{'Control010'}, {'Control011'}, {'Control089'},{'Control046'}, {'Control121'}, {'Control141'}, {'Control186'}];

%Set the output file location
outLocation = 'C:\Users\Tim\Desktop\Nif Project Figures\Motif Continuity Analysis\';
outName = 'Lesion Summary MC 051815.mat';

%Ask user for data directory
dataFolder = uigetdir('C:\Users\Tim\Desktop\Nif Project Figures\Motif Continuity Analysis\');

%Get all data files that match the criterion
crit = '*MCanal.mat';
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
%     parsed(i).files = dataOut.fileSets;

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
    parsed(i).dayRaw = mcSet(:,1);
    
    %Separate out day and night measurements
%     rems = dataOut.time-floor(dataOut.time);
%     mornIdx = rems <= (14/24);
%     nightIdx = rems > (14/24);
%     t = floor(dataOut.time); %Round down to nearest 0.5
%     t(nightIdx) = t(nightIdx) + 0.5;

    for j = 1:size(mcSet,1)
        %Unpack the time to day/night
        parsed(i).dayInt(2*j-1) = mcSet(j,1);
        parsed(i).dayInt(2*j) = mcSet(j,1)+0.5;
        
        %Unpack MC to day/night
        parsed(i).MC(2*j-1) = mcSet(j,2);
        parsed(i).MC(2*j) = mcSet(j,3);
    end
    
    ind = ~isnan(parsed(i).MC);
    parsed(i).dayInt = parsed(i).dayInt(ind);
    parsed(i).MC = parsed(i).MC(ind);
    
    if strcmp(parsed(i).type, 'control')
        parsed(i).dayInt = parsed(i).dayInt + 7; %all controls come before the lesion day
    else
%         parsed(i).dayInt = parsed(i).dayInt(parsed(i).dayInt~=0.5);
%         parsed(i).MC = parsed(i).MC(parsed(i).dayInt~=0.5);
    end
    
    %Clear the loaded data
    clear('mcSet')
    
end


%Organize the parsed data for plotting
xs = 1:9;
plotData = NaN(length(parsed), length(xs));
for i = 1:length(parsed)
    %Capture the baseline lengths
    indx = find(parsed(i).dayInt <= 0);
    meanPre = nanmean(parsed(i).MC(indx));
    
    plotData(i,1) = meanPre;
    
    %Sort the remaining data days
    critPnts = [0.5, 1, 2, 3, 4, 5, 6, 7];
    for j = 1:length(critPnts)
        if ismember(critPnts(j),parsed(i).dayInt)
            indx = find(parsed(i).dayInt == critPnts(j),'1', 'first');
            plotData(i,j+1) = parsed(i).MC(indx);
%         elseif j ~= 1 && ismember(critPnts(j)+.5,parsed(i).dayInt) %if there is no morning data, take the night if present
        elseif ismember(critPnts(j)+.5,parsed(i).dayInt) %if there is no morning data, take the night if present
            indx = find(parsed(i).dayInt == critPnts(j)+.5,'1', 'first');
            plotData(i,j+1) = parsed(i).MC(indx);
        end
    end

    plotDataNorm(i,:) = plotData(i,:)./plotData(i,1);

end

%Parse the data for day/night recovery
% isUni = strcmp(type,'unilateral');
% uniSub = parsed(isUni);
% for i = 1:length(uniSub)
%     %Separate into pre-post lesion and day/night
%     preIndx = find(uniSub(i).dayInt <= 0);
%     dayIndx = find(uniSub(i).dayInt == round(uniSub(i).dayInt));
%     
%     
%     meanPre = (parsed(i).MC(indx));
% 
% end


%Retrieve name and lesion type fromt he parsing structure
name = getStructField(parsed,'name');
type = getStructField(parsed,'type');

%Create plotting indices
isBi = strcmp(type,'bilateral');
isUni = strcmp(type,'unilateral');
isC = strcmp(type,'control');

figure(1); clf
offset = 0.0125;
p = plotDataNorm(isBi,:);
biMean = nanmean(p,1); biSem = nanstd(p,1)./sqrt(size(p,1));
bxs = [1,4:9]; %bilateral subset
plot(bxs-offset,nanmean(p(:,bxs),1),'.b'); hold on
errorbar(bxs-offset,nanmean(p(:,bxs),1), nanstd(p(:,bxs),1)./sqrt(size(p,1)),'-b')

p = plotDataNorm(isUni,:);
uniMean = nanmean(p,1); uniSem = nanstd(p,1)./sqrt(size(p,1));
uxs = [1:9]; %unilateral subset
plot(uxs+offset,nanmean(p(:,uxs),1),'.r')
errorbar(uxs+offset,nanmean(p(:,uxs),1), nanstd(p(:,uxs),1)./sqrt(size(p,1)),'-r')

p = plotDataNorm(isC,:);
cMean = nanmean(p,1); cSem = nanstd(p,1)./sqrt(size(p,1));
cxs = [1,3:9];
plot(cxs,nanmean(p(:,cxs),1),'.k')
errorbar(cxs,nanmean(p(:,cxs),1), nanstd(p(:,cxs),1)./sqrt(size(p,1)),'-k')
xlim([0.5, 9.5]); ylim([0.5,1.1]);
xlabel('Time Post-Lesion', 'FontSize', 10); ylabel('Normalized Motif Truncation (%)', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10, 'XTick', 1:9, 'XTickLabels', [{'Pre'},{'1hr'},{'1d'}, {'2d'}, {'3d'}, {'4d'}, {'5d'}, {'6d'}, {'7d'}], 'YTick', [0.5:0.25:1])

set(gcf, 'Units', 'Inches','Position', [0,0,4,4])


[t.PreUni, p.PreUni, ~] = ttest(plotDataNorm(isUni,1),plotDataNorm(isUni,4));      %p = 0.1819
[t.PreBi, p.PreBi, ~] = ttest(plotDataNorm(isBi,1),plotDataNorm(isBi,5));            %p = 0.1136

[t.UniPreH1, p.UniPreH1, ~] = ttest(plotDataNorm(isUni,1),plotDataNorm(isUni,2));      %p = 0.039
[t.BiPreD2, p.BiPreD2, ~] = ttest(plotDataNorm(isBi,1),plotDataNorm(isBi,4))     %p = 0.1245


%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model for unilateral
q = plotDataNorm(isUni,:);
t = table(q(:,1), q(:,2), q(:,3), q(:,4), q(:,5), q(:,6), q(:,7), q(:,8), q(:,9),'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9'});
Meas = table([1 2 3 4 5 6 7 8 9]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas9~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm);
epsilon(rm)

%%%%%%%
%Results for unilateral 
%%%%%%%
% ranovatbl = 
%                                                   SumSq     DF     MeanSq        F         pValue         pValueGG   pValueHF   pValueLB
%                                                   _______    __    _________    ______    __________    ________    ________    ________
% 
%     (Intercept):Measurements         0.44421     8     0.055527    7.0133    0.00050279    0.090782    0.029485    0.1179  
%     Error(Measurements)               0.12668    16    0.0079173

%     W    ChiStat    DF    pValue
%     _    _______    __    ______
%     0    -Inf       35    1 

%Do multiple comparisons
%multcompare(rm, 'Measurements', 'ComparisonType', 'tukey-kramer')

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
%     Pre                                 1h                5.044e-04  
%     Pre                                 1d             	0.1165
%     Pre                                 2d                0.9998
%     Pre                                 3d                0.9988
%     Pre                                 4d                1.0000  
%     Pre                                 5d                1.0000  
%     Pre                                 6d                1.0000  
%     Pre                                 7d                0.9992  

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
%                                                SumSq      DF     MeanSq        F        pValue      pValueGG    pValueHF   pValueLB
%                                             ________    __    _________    ______    _______    ________    _______    ________
%     (Intercept):Measurements    0.057658     6    0.0096097    4.0724    0.0093823    0.11015     0.071706    0.13692 
%     Error(Measurements)          0.042474    18    0.0023597    
%
%       W    ChiStat    DF    pValue
%       _    _______    __    ______
% 
%       0    Inf        20    0


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
%     Pre                                 2d                0.0088
%     Pre                                 3d                0.8327
%     Pre                                 4d                0.9607  
%     Pre                                 5d                0.9967  
%     Pre                                 6d                0.9192  
%     Pre                                 7d                0.9468  


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

%save([outLocation, outName], 'output')





