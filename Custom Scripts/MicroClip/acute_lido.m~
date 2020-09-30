function acute_lido
%This function will take as input a user-specified source folder and range
%of files and do some plotting and summary stats on the selected files.
%
%
% Updated by TMO 9/29/20
clear

%Request input from the user on which folder to process
folder = '/Users/tim/Desktop/Acute Recordings/LLR20';

%A few current sweeps from LLR20
saline_start = 150448;
saline_end = 150505;

lido_start = 150648;
lido_end = 150652;

washout_start = 151746;
washout_end = 151813;

%Processing constants
respMin = 0.75;
respMax = 6;
select = [1, 2, 5, 6];
sym = {'o', 'x', 'sq', 'd'};

%Get the folder file contents and check for empties
files = dir([folder, filesep, 'stim*.mat']);
if isempty(files)
    disp('Nothing in the folder, bro... try again.')
    return
end

%Extract timestamp from the filenames
tStamp = [];
for i = 1:numel(files)
    s = files(i).name;
    sp = regexp(s, '_', 'split');
    sps = sp{end}(1:6);
    tStamp(i) = str2double(sps);
end

%Derive selected files from masking on time bounds
saline_mask = (tStamp >= saline_start) & (tStamp <= saline_end);
lido_mask = (tStamp >= lido_start) & (tStamp <= lido_end);
washout_mask = (tStamp >= washout_start) & (tStamp <= washout_end);

%Data collection
[saline_responses, times] = collect_data(folder, files(saline_mask));
lido_responses = collect_data(folder, files(lido_mask));
washout_responses = collect_data(folder, files(washout_mask));

%Process data
saline_CMS = process_data(saline_responses);
lido_CMS = process_data(lido_responses);
washout_CMS = process_data(washout_responses);

%Plot the data
f1 = figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 10, 8])
plot_data(f1, times, saline_CMS, select, 'k')
plot_data(f1, times, lido_CMS, select, 'r')
plot_data(f1, times, washout_CMS, select, 'g')

%Calculate Vpp
mask = times >= respMin & times <= respMax;
saline_vpp = vpp_data(saline_CMS, mask);
lido_vpp = vpp_data(lido_CMS, mask);
washout_vpp = vpp_data(washout_CMS, mask);

%Plot the vpp across conditions
chunk = [saline_vpp, lido_vpp, washout_vpp];
mean_vpp = mean(chunk(select, :), 1);
std_vpp = std(chunk(select, :), 1);

f2 = figure(2); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 10, 5, 6])
bar(1:3, mean_vpp); hold on
errorbar(1:3, mean_vpp, std_vpp, '.k')
for i = 1:numel(select)
    jit = (rand(1,3)-0.5)/20;
    plot([1:3] + jit, chunk(select(i), :), [sym{i}, 'k'])
end
xlim([0.5, 3.5])
ylim([0, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.2:0.4, 'XTickLabel', {'Saline', 'Lidocaine', 'Washout'})
ylabel('Vpp (mV)')

% Compute significance
dataMat = chunk(select,:);

%Check that all data are from normal distributions (for t-tests)
h_ksStats = []; p_ksStats = [];
for i = 1:3
        [h_ksStats(i), p_ksStats(i)] = kstest(dataMat(:,i));
end
allNormal = all(h_ksStats(:)); % allNormal = false ==> lidocaine is not normal.


%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model
t = table(dataMat(:,1), dataMat(:,2), dataMat(:,3),'VariableNames',{'meas1','meas2','meas3'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1', 'WithinDesign', Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

%Do multiple comparisons
mc = multcompare(rm, 'Measurements', 'ComparisonType', 'tukey-kramer')

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
ma = mauchly(rm)
ep = epsilon(rm)

%Do Dunnett's Test with PBS as the control
f = table2array(ranovatbl);
stats = [];
stats.means = mean(dataMat,1);
stats.df = f(2,2);
stats.n = [3,3,3];
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, [2, 3], 1 )


% ranovatbl =
% 
%                                SumSq     DF     MeanSq        F        pValue      pValueGG    pValueHF    pValueLB
%                                 _______    __    _________    ______    _________    ________    ________    ________
% 
%     (Intercept):Measurements    0.09484    2       0.04742    12.651    0.0070429    0.035453    0.031977    0.037906
%     Error(Measurements)         0.02249    6     0.0037484                                                           
%
%
% mc = 
%     Measurements_1    Measurements_2    Difference     StdErr      pValue       Lower        Upper  
%     ______________    ______________    __________    ________    ________    _________    _________
% 
%           1                 2             0.21465     0.057735    0.067335    -0.026613      0.45591
%           1                 3            0.075541      0.01506    0.030819     0.012611      0.13847
%           2                 1            -0.21465     0.057735    0.067335     -0.45591     0.026613
%           2                 3            -0.13911     0.045413     0.10746     -0.32888     0.050665
%           3                 1           -0.075541      0.01506    0.030819     -0.13847    -0.012611
%           3                 2             0.13911     0.045413     0.10746    -0.050665      0.32888
% 
%
% ma = 
%        W        ChiStat    DF     pValue 
%     ________    _______    __    ________
% 
%     0.074078    5.2053     2     0.074078
%   
%
%     Uncorrected    GreenhouseGeisser    HuynhFeldt    LowerBound
%     ___________    _________________    __________    __________
% 
%          1              0.51923          0.54902         0.5    
%          
% pvals =
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
%  
%     PBS                   Lido             0.0092 
%     PBS                   Wash             0.2941



%%%%%%%%%%%%%
%Subfunctions

function vpp = vpp_data(responses, mask)
vpp = [];
for i = 1:6
    mean_trace = mean(responses(:,:,i))*10e3;
    snip = mean_trace(mask);
    vpp = [vpp; max(snip) - min(snip)];
end


function CMS_responses = process_data(responses)

%Demean
for i = 1:6
    %find offset
    delta = (mean(responses(:,1:75,i),2));

    responses(:,:,i) = responses(:,:,i) - delta;
end

%CMS
for i = 1:6
    % Define averaging mask
    mask = true(1, 6);
    mask(i) = false;

    %Process trials separately
    for j = 1:size(responses,1)
        %Define common mode
        commonMode = mean(responses(j, :, mask), 3);
        
        %Subtract common from signal
        CMS_responses(j,:,i) = responses(j,:,i) - commonMode;
    end
end


function plot_data(f, times, responses, select, col)
figure(f)

%Blank artifact
T = times;
preT = find(times < -0.05);
postT = find(times > 0.75);

for j = 1:numel(select)
    i = select(j);
    times = T;
    %Plot the mean trace
    subplot(2, 2, j)
    shadedErrorBar(times(preT), mean(responses(:,preT,i))*10e3, std((responses(:,preT,i))*10e3, 1, 1), col, 1); hold on
    if strcmp(col, 'g')
        times = times-0.2;
        postT = find(times > 0.75);
    end
    shadedErrorBar(times(postT), mean(responses(:,postT,i))*10e3, std((responses(:,postT,i))*10e3, 1, 1), col, 1); hold on
    
    %Format the figure
    axis tight
    xlim([-2, 6])
    ylim([-0.35, 0.35])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -2:2:6, 'YTick', -0.35:0.35:0.35)
    ylabel('Evoked Response (mV)')
    xlabel('Time (ms)')
    
end


function [responses, times] = collect_data(folder, files)
responses = [];
for i = 1:numel(files)
    cur_file = files(i).name;
    load([folder, filesep, cur_file]);
    
    responses = cat(1, responses, data.tdt.response);
end
times = data.tdt.times_aligned*1000;




















