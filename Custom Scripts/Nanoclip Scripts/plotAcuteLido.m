%Plot and calculate significance of the Vpp measurements from the acute
%nanoclip recordings. 
%
%Written by TMO 05/24/19

%% Clear the workspace
clear all

%Source and destination data
sLoc = '/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/Acute Recording PP.mat';

%Load variables from file
load(sLoc);
if ~exist('acuteRecData', 'var')
    disp('Uh-oh... something is wrong')
    return
end

%Normalize all measurements to baseline PP
normData = [];
numBirds = numel(acuteRecData);
for i = 1:numBirds
    base(i) = acuteRecData(i).baselinePP./acuteRecData(i).baselinePP;
    saline(i) = acuteRecData(i).salinePP./acuteRecData(i).baselinePP;
    lido(i) = acuteRecData(i).lidocainePP./acuteRecData(i).baselinePP;
    washout(i) = acuteRecData(i).washoutPP./acuteRecData(i).baselinePP;
end

% %Format for easy stats/plotting
% base = getFieldVector(normData, 'baseline');
% saline = getFieldVector(normData, 'saline');
% lido = getFieldVector(normData, 'lidocaine');
% washout = getFieldVector(normData, 'washout');

dataMat = [base', saline', lido', washout'];

%Descriptive stats
mPP = mean(dataMat,1);
sPP = std(dataMat,1,1);

%% Plot the bar plot 
figure(41); clf
sym = {'o', 'sq', '*'};
cats = 1:numel(mPP);
labels = {'Baseline'; 'Saline'; 'Lidocaine'; 'Washout'};
bar(cats, mPP, 0.6, 'LineStyle', 'none'); hold on
errorbar(cats, mPP, sPP, 'LineStyle', 'none', 'Marker', '.', 'Color', 'k', 'LineWidth', 1.5)
for i = 1:3
    scatter([1, 2, 3, 4], [base(i), saline(i), lido(i), washout(i)], 50, 'k', sym{i})
end

% axis square
xlim([0.5, 4.5]); ylim([0, 1.05]);
ylabel('Normalized V_p_p')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:4, 'XTickLabels', labels, 'YTick', [0,1])
set(gcf, 'Units', 'Inches', 'Position', [10.5, 5.5, 6, 5.5])





%% Compute significance tests

%Check that all data are from normal distributions (for t-tests)
h_ksStats = []; p_ksStats = [];
for i = 1:4
        [h_ksStats(i), p_ksStats(i)] = kstest(dataMat(:,i));
end
allNormal = all(h_ksStats(:)); % allNormal = false ==> lidocaine is not normal.

% %Failed the normality test, so performing Wilcox rank sum test
% h_wrsStats = []; p_wrsStats = [];
% for i = 1:3
%         [p_wrsStats(i), h_wrsStats(i)] = ranksum(dataMat(:,2),dataMat(:,i+1));
% end

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model
t = table(dataMat(:,1), dataMat(:,2), dataMat(:,3), dataMat(:,4),'VariableNames',{'meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas4~1', 'WithinDesign', Meas);

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
stats.n = [3,3,3,3];
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, [1, 3, 4],2 )


% ranovatbl =
%                                  SumSq      DF     MeanSq        F         pValue      pValueGG    pValueHF     pValueLB 
%                                 ________    __    _________    ______    __________    ________    _________    _________
% 
%     (Intercept):Measurements      1.7789    3       0.59298    189.02    2.5007e-06    0.004736    0.0034444    0.0052487
%     Error(Measurements)         0.018823    6     0.0031371                                                              
%
%
%
% mc = 
%     Measurements_1    Measurements_2    Difference     StdErr       pValue       Lower       Upper  
%     ______________    ______________    __________    ________    __________    ________    ________
% 
%           1                 2            0.024479     0.028356       0.82626    -0.17197     0.22093
%           1                 3             0.93292     0.011478    0.00037882      0.8534      1.0124
%           1                 4             0.12737     0.071558       0.46846     -0.3684     0.62313
%           2                 1           -0.024479     0.028356       0.82626    -0.22093     0.17197
%           2                 3             0.90844     0.021086     0.0013472     0.76235      1.0545
%           2                 4             0.10289     0.043542       0.32375    -0.19878     0.40455
%           3                 1            -0.93292     0.011478    0.00037882     -1.0124     -0.8534
%           3                 2            -0.90844     0.021086     0.0013472     -1.0545    -0.76235
%           3                 4            -0.80555     0.064432      0.015817     -1.2519    -0.35916
%           4                 1            -0.12737     0.071558       0.46846    -0.62313      0.3684
%           4                 2            -0.10289     0.043542       0.32375    -0.40455     0.19878
%           4                 3             0.80555     0.064432      0.015817     0.35916      1.2519
%
%
% ma = 
%     W    ChiStat    DF    pValue
%     _    _______    __    ______
% 
%     0      Inf      5       0   
%   
%
% ep = 
%     Uncorrected    GreenhouseGeisser    HuynhFeldt    LowerBound
%     ___________    _________________    __________    __________
% 
%          1              0.34202          0.36901       0.33333  
%          
% pvals =
% 
%     0.9077    0.0000    0.1474
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
% 
%     PBS                   Base             0.9077    
%     PBS                   Lido             6.39e-06 
%     PBS                   Wash             0.1474








