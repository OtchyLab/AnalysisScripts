function acute_lido_summary
%This function will take as input the output of function acute_lido.m, and
%create a summary plot and across-animal statistics. Figures to be used in
%Rowan et al (2020) paper.
%
%
% Updated by TMO 10/05/20

%Load the summary stats from file
folder = '/Users/tim/Desktop/Acute Recordings';
file = 'Lidocaine Exp Summary Data.mat';
load([folder, filesep, file])

%Plot the vpp across conditions
chunk = [saline_vpp, lido_vpp, washout_vpp];
mean_vpp = mean(chunk, 1);
std_vpp = std(chunk, 1);

f = figure(1); clf
sym = {'o', 'x', 'sq', 'd'};
set(gcf, 'Units', 'Inches', 'Position', [10, 10, 5, 6])
bar(1:3, mean_vpp); hold on
errorbar(1:3, mean_vpp, std_vpp, '.k')
for i = 1:size(chunk,1)
    jit = (rand(1,3)-0.5)/20;
    plot([1:3] + jit, chunk(i, :), [sym{i}, 'k'])
end
xlim([0.5, 3.5])
ylim([0, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.2:0.4, 'XTickLabel', {'Saline', 'Lidocaine', 'Washout'})
ylabel('Vpp (mV)')

% Compute significance
dataMat = chunk;

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
