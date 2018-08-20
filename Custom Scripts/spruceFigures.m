% Script is to spruce up final figures for the Lesions/Inactivations paper
% TMO 04-28-2015

%% Figure 2 Pharma Lesion Feature Heatmaps
%To execute this, load the PharmaLesion Summary Stats file that contains the pharmaStats structure. The example in the figure
%uses Grn011, which is entry #2 in the structure. This cell code will extract the heatmaps for the requested days, up-sample
%them, trim to size, plot and format.

%Get images/pdfs
pre = squeeze((pharmaStats(2).Nsm(11,:,:)));
post =  squeeze((pharmaStats(2).Nsm(9,:,:)));

%Resize and smooth with bicubic interpolation
% preBig = imresize(pre, 2);
% postBig = imresize(post,2);
preBig = pre;
postBig = post;

%Limit the range to 0-300ms from the original 0-350ms
[r, c] = size(preBig);
subRows = floor(r*(300/350)); %reduce timescale to 300/350

preBig_tr = preBig(1:subRows,:);
postBig_tr = postBig(1:subRows,:);

ratio = (subRows/300);

%Plot to figure
load('inactCmap.mat');
figure(3); clf
subplot(2,1,1); cla
imagesc(preBig_tr', [0, 0.0018]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
title('Grn011 Pre-Post Heatmaps')

subplot(2,1,2); cla
imagesc(postBig_tr', [0, 0.0018]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([1,100,200,300]), 'XTickLabels', [0,100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

set(gcf, 'Units', 'Inches', 'Position', [0,0,3,4])


%% Figure 2 Lesion Distances Plot
%To execute this, load the Summary Pharma Data (produced by pharmaSumStatsPlot.m) file that contains all the data from
%that script (not just the minimal output structure.
%This cell code will plot the distances between heatmaps in all conditions and format it to be similar to Julie's plots

%Plot to figure
% Manually alter scaling of the jet colormap so that the 2nd blue tick is aligned with the 2nd color bin
figure(2); clf
q = [];
% markers = [{'o'}, {'s'}, {'d'}, {'^'},{'p'},{'h'}];

b = bar(1:4, mean(subdist,1), 0.3, 'FaceColor', 'none'); hold on
for i = 1:5
    q(i) = plot(1:4, subdist(i,:), '.k');
    set(q(i), 'Marker', 'o', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'LineStyle', 'none');
end
eb = errorbar(1:4, mean(subdist,1), std(subdist,1)./sqrt(size(subdist,1)),'.k');
xlim([0.5,4.5]); ylim([0,50]);
ylabel('Wasserstein Distance', 'FontSize', 12)
title(['Summary Lesion'], 'FontSize', 12);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4], 'XTickLabel', [{'Control'}, {'Day 2'}, {'Day 3'}, {'Day 8'}], 'YTick', [0:25:75]);


set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5])

[t.C2, p.C2, ~] = ttest2(subdist(:,1), subdist(:,2) ); % p = 0.1981
[t.C3, p.C3, ~] = ttest2(subdist(:,1), subdist(:,3) ); % p = 0.4693
[t.C7, p.C7, ~] = ttest2(subdist(:,1), subdist(:,4) );  % p = 0.8568

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model
t = table(subdist(:,1), subdist(:,2), subdist(:,3), subdist(:,4),'VariableNames',{'meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas4~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

%Do multiple comparisons
multcompare(rm, 'Measurements', 'ComparisonType', 'tukey-kramer')

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm)

%   ranovatbl = 
%                                               SumSq     DF    MeanSq      F         pValue   pValueGG   pValueHF    pValueLB
%                                               ______     __    ______    ______    _______   ________    ________    ________
% 
%     (Intercept):Measurements      17.055     3       5.6849    2.2065   0.14008    0.1717         0.14008       0.21161 
%     Error(Measurements)             30.917    12      2.5764          
%
%Note that these values are identical to the output of the anova_rm.m script (from MFE) which is a simple rmanova, suggesting
%that the setup for the modeling is correct.
%        
% BECAUSE DOES NOT PASS ANOVA, NO NEED FOR POST-HOC
%
%        W             ChiStat    DF    pValue 
%     ________    _______    __    _______
%     0.058632      7.7215      5     0.1722
    
%% Figure 2 Inactivation Feature Heatmaps
%To execute this, load the Summary Inactivation Data (produced by inactSumStatsPlotV2.m) file that contains all the data from
%that script (not just the minimal output structure. The example in the figure uses Grn115, which is entry #4 in the input structure.
%This cell code will recalculate the  heatmaps for the conditions, up-sample them, trim to size, plot and format.

%This selects which bird to process
i = 4; 
birdMask = strcmp(birds(i), name);

%Generates all the mean heatmaps for all conditons
preNmean = squeeze(mean(preN(birdMask, :, :),1));
postNmean = squeeze(mean(postN(birdMask, :, :),1));
inactNmean = squeeze(mean(injN(birdMask & inactMask, :, :),1));
PBSNmean = squeeze(mean(injN(birdMask & PBSMask, :, :),1));
elevNmean = squeeze(mean(injN(birdMask & elevMask, :, :),1));

%Rescale all for hi-res images
% preBig = imresize(preNmean, 10);
% inactBig = imresize(inactNmean, 10);
% PBSBig = imresize(PBSNmean, 10);
% elevBig = imresize(elevNmean, 10);
% postBig = imresize(postNmean, 10);
preBig = preNmean;
inactBig = inactNmean;
PBSBig = PBSNmean;
elevBig = elevNmean;
postBig = postNmean;

%Trim to 300ms from 350ms
[r, c] = size(preBig);
subRows = floor(r*(300/350)); %reduce timescale to 300/350
ratio = (subRows/300);

preBig_tr = preBig(1:subRows,:);
inactBig_tr = inactBig(1:subRows,:);
PBSBig_tr = PBSBig(1:subRows,:);
elevBig_tr = elevBig(1:subRows,:);
postBig_tr = postBig(1:subRows,:);

%Plot to figure
load('inactCmap.mat');
figure(5); clf
subplot(5,1,1); cla
imagesc(preBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
title('Grn115 Heatmaps')

subplot(5,1,2); cla
imagesc(inactBig_tr', [0, 0.0028]); axis xy; %colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(5,1,3); cla
imagesc(PBSBig_tr', [0, 0.0028]); axis xy; %colormap(heatCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(5,1,4); cla
imagesc(elevBig_tr', [0, 0.0028]); axis xy; %colormap(heatCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(5,1,5); cla
imagesc(postBig_tr', [0, 0.0028]); axis xy; %colormap(heatCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([1,100,200,300]), 'XTickLabels', [0,100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
xlabel('Syl Dur (ms)'); ylabel('log(Entropy)')

set(gcf, 'Units', 'Inches', 'Position', [0,0,3,8])

%% Figure 2Inactivation Distances Plot
%To execute this, load the Summary Inactivation Data (produced by inactSumStatsPlotV2.m -- right now, named
%"Summary Inactivation Data Final 05162015.mat") file that contains all the data from
%that script (not just the minimal output structure.
%This cell code will plot the distances between heatmaps in all conditions and format it to be similar to Julie's plots

%Plot to figure
% Manually alter scaling of the jet colormap so that the 2nd blue tick is aligned with the 2nd color bin
figure(2); clf
q = [];
%markers = [{'o'}, {'s'}, {'d'}, {'^'},{'p'},{'h'}];
b = bar(1:5, [mean(ddDist), mean(inactDist), mean(PBSDist), mean(elevDist), mean(ppDist)], 0.3, 'FaceColor', 'none'); hold on
for i = 1:5
    q(i) = plot(1:5, [ddDist(i), inactDist(i), PBSDist(i), elevDist(i),ppDist(i)], '.k');
    set(q(i), 'Marker', 'o', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'LineStyle', 'none');
end
eb = errorbar(1:5, [mean(ddDist), mean(inactDist), mean(PBSDist), mean(elevDist), mean(ppDist)], [std(ddDist), std(inactDist), std(PBSDist), std(elevDist), std(ppDist)]./sqrt(5),'.k');
xlim([0.5,5.5]); ylim([0,75]); %axis square
ylabel('Wasserstein Distance', 'FontSize', 12)
title(['Summary Inactivation'], 'FontSize', 12);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4, 5], 'XTickLabel', [{'Control'}, {'Muscimol'}, {'Saline'}, {'Elevated'}, {'Washout'}], 'YTick', [0:25:75], 'LineWidth', 2, 'FontSize', 10)   

set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5])

%Hypothesis Testing
[testOut.ddpp,pOut.ddpp,ci.ddpp,~] = ttest([ddDist],[ppDist]);                               %p = 0.6423
[testOut.ddInact,pOut.ddInact,ci.ddInact,~] = ttest([ddDist],[inactDist]);                 %p = 0.0067
[testOut.ddPBS,pOut.ddPBS,ci.ddPBS,~] = ttest([ddDist],[PBSDist]);                  %p = 0.6532
[testOut.ddElev,pOut.ddElev,ci.ddElev,~] = ttest([ddDist],[elevDist]);                      %p = 0.5261
[testOut.inactPBS,pOut.inactPBS,ci.inactPBS,~] = ttest([inactDist],[PBSDist]);    %p = 0.0084
[testOut.inactElev,pOut.inactElev,ci.inactElev,~] = ttest([inactDist],[elevDist]);        %p = 0.0160
[testOut.PBSelev,pOut.PBSelev,ci.PBSelev,~] = ttest([PBSDist],[elevDist])            %p = 0.0378

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model
t = table(ddDist, inactDist, PBSDist, elevDist,ppDist,'VariableNames',{'cond1', 'cond2', 'cond3', 'cond4', 'cond5'});
Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
rm = fitrm(t,'cond1-cond5~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm)

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm)
epsilon(rm)
% ranovatbl = 
%                                              SumSq     DF    MeanSq      F         pValue      pValueGG       pValueHF      pValueLB
%                                               ______   __    ______    ______    __________    _________    _________    ________
% 
%     (Intercept):Measurements    4153.7     4    1038.4    19.722    4.9997e-06    0.0071906    0.0043469    0.011325
%     Error(Measurements)          842.43    16    52.652    
% 
%         W               ChiStat    DF     pValue 
%     __________    _______    __    ________
%     0.00018689       20.747     9     0.013822
% 
%     Uncorrected    GreenhouseGeisser    HuynhFeldt    LowerBound
%     ___________    _________________    __________    __________
% 
%             1                     0.29219              0.33942            0.25      


%Do Dunnett's Test with PBS as the control
f = table2array(ranovatbl);
stats = [];
stats.means = mean([ddDist, inactDist, PBSDist, elevDist,ppDist],1);
stats.df = f(2,2);
stats.n = [size(ddDist,1), size(inactDist,1), size(PBSDist,1), size(elevDist,1), size(ppDist,1),];
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, [1, 2, 4, 5],3)

% pvals = 
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
% 
%     PBS                                 Base             0.9984    
%     PBS                                 Musc             1.013e-05  
%     PBS                                 Elev               0.9237
%     PBS                                 Wash            1.0000   


%% Figure 2 Inactivation and HVC lesion PDFs
%To execute this, load the Summary Inactivation Data (produced by inactSumStatsPlotV2.m) file that contains all the data from
%that script (not just the minimal output structure. Also need to load into memory the duration lists for the five birds from
%Aronov

%This cell code will plot mean traces for the Nif inactivation and HVC lesion experiments

%Capture the mean 1-D pdfs for the Nif inactivations
inactY = [];
for i = 1:length(birds)
    %Masks
    birdMask = strcmp(birds(i), name);
    inactMask = strcmp(conditions(1), type);

    %Get the mean trace
    inactY(i,:) = mean(injY(birdMask & inactMask, :),1);
end

%Names from Aronov are lesionDur129, 150, 165, 189, 213
binSize = 5; minval = 10; maxval = 400;
lesionY = [];
[~, lesionY(1,:)] = epdf_cbins(lesionDur129,binSize,minval,maxval);
[~, lesionY(2,:)] = epdf_cbins(lesionDur150,binSize,minval,maxval);
[~, lesionY(3,:)] = epdf_cbins(lesionDur165,binSize,minval,maxval);
[~, lesionY(4,:)] = epdf_cbins(lesionDur189,binSize,minval,maxval);
[~, lesionY(5,:)] = epdf_cbins(lesionDur213,binSize,minval,maxval);

%Calculate the mean trace for each dataset
inactYmean = mean(inactY,1);
lesionYmean = mean(lesionY,1);
lesionYmean(1) = inactYmean(1);

%Smooth it with a sliding window
inactSmooth = nanmean(windowTS(inactYmean, 5, 1, 'pad', 'boxcar')');
lesionSmooth = nanmean(windowTS(lesionYmean, 5, 1, 'pad', 'boxcar')');

%Plot only the mean of each group
figure(52); clf
plot(10:5:400, mean(inactSmooth,1),'r', 'LineWidth', 2, 'DisplayName', 'Nif Inactivation (n = 5)'); hold on
plot(10:5:400, mean(lesionSmooth,1), 'b', 'LineWidth', 2, 'DisplayName', 'HVC Lesion (n = 5)');
xlim([0, 300]); ylim([0, 0.07])
xlabel('Syllable Duration (ms)'); ylabel('P(t)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [0, 100, 200, 300], 'YTick', [0, 0.05], 'LineWidth', 2, 'FontSize', 10)
legend('show'); legend('boxoff')
set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5])

%Similarity meansures:
inactLesionCorr = corr(mean(inactSmooth,1)',mean(lesionSmooth,1)') % on 4/29/2015, correlation was 0.9420

[h, p] = kstest2(mean(inactSmooth,1),mean(lesionSmooth,1)) %p = 0.4007 -- not significantly different

%% Figure 1 Proper Statistics for Julie's Inactivation Data
%This script operates on the data Julie emailed to TMO on 8/13/2015. Those numbers are hardcoded below. The output of this
%script are simily the ANOVA and Dunnett's Post-Hoc Test Results. (She handled the plotting of the shown figure elsewhere.)

%Hardcoded data from Julie; ordering [baseline, PBS, musc]
Kansas =  [0.8205, 0.8629, 0.1127];
Iraq = [0.6416, 0.5629, 0.2656];
Anax = [0.5745, 0.6173, 0.1445];
Iran = [0.6339, 0.5481, 0.2391];
Jane = [0.6905, 0.7657, 0.0833];

%Concatenate into single block
test = [Kansas; Iraq; Anax; Iran; Jane];

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/13/15)
%Establish repeated measures model
t = table(test(:,1), test(:,2), test(:,3),'VariableNames',{'cond1', 'cond2', 'cond3'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'cond1-cond3~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm)

%Test the equality of variances
vartestn(test)

% Test sphericity assumption using mauchly test -- if p<0.05, there are siginificant differences in variance between groups and a correction must be applied
mauchly(rm)
epsilon(rm)

% ranovatbl = 
%                                               SumSq      DF    MeanSq       F          pValue      pValueGG     pValueHF     pValueLB 
%                                             ________    __    _______    ______    _______     _________    _________    _________
% 
%     (Intercept):Measurements     0.84253    2     0.42126    35.731    0.00010273    0.0030133    0.0022792    0.0039362
%     Error(Measurements)         0.094318    8     0.01179          
%
%        W       ChiStat    DF     pValue 
%     _______    _______    __    ________
% 
%     0.13423    6.0246     2     0.049179
% 
%     Uncorrected    GreenhouseGeisser    HuynhFeldt    LowerBound
%     ___________    _________________    __________    __________
% 
%             1                     0.53597                 0.57371                 0.5


%Do Dunnett's Test with PBS as the control
f = table2array(ranovatbl);
stats = [];
stats.means = mean(test,1);
stats.df = f(2,2);
stats.n = [size(test,1), size(test,1), size(test,1)];
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, [1, 3], 2)

% pvals = 
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
% 
%     PBS                                 Base              0.9999
%     PBS                                 Musc             0.0002

%% Figure 3 Steffen's Heatmaps and Distance Plot
%To execute this, load all the data from Steffen's Analysis Summary (named "All Analysis Data 0422.mat" for now, and produced
%by RatAnalysisSW.m) file that contains all the data from that script (not just the minimal output structure.
%
%This cell script Plots and example heatmap and a summary distance figure

%Select the data for the example rat. I'm choosing Rat #2 (no idea his real "name"
i = 2; 
onHeat = squeeze(onPDF(i,:,:));
offHeat = squeeze(mean(squeeze(offPDF(i,:,:,:)),1));

%Resize
onHeatBig = imresize(onHeat, 10);
offHeatBig = imresize(offHeat, 10);
numBins = size(onHeatBig,1);

%Plot Heatmaps
% Manually alter colomap by dragging secnd blue tick to Index = 2, and the third to Index = 11
% figure(9); clf
% subplot(2,1,1)
% imagesc(offHeatBig'); cs = caxis; colormap(jet); axis xy; axis image
% set(gca,'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'XTick', linspace(1,numBins,3), 'YTick', linspace(1,numBins,3),'XTickLabel', 3:2:8, 'YTickLabel',3:2:8)
% 
% subplot(2,1,2)
% imagesc(onHeatBig', cs); axis xy; axis image
% set(gca,'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'XTick', linspace(1,numBins,3), 'YTick', linspace(1,numBins,3),'XTickLabel', 3:2:8, 'YTickLabel',3:2:8)
% 
% xlabel('log(Tap1 Dur)'); ylabel('log(IPI-Tap1 Dur)')
% set(gcf,'Units','Inches','Position',[0,0,3,4])

%Plot the summary distances
figure(10); clf
q = [];
markers = [{'o'}, {'s'}, {'d'}, {'^'},{'p'},{'h'}];
b = bar(1:2, [mean(offoffEMDm), mean(onoffEMDm)], 0.3, 'FaceColor', 'none'); hold on
for i = 1:2
    q(i) = plot(1:2, [offoffEMDm(i), onoffEMDm(i)],'.k');
    set(q(i), 'Marker', markers{i}, 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'LineStyle', 'none');
end
eb = errorbar(1:2, [mean(offoffEMDm), mean(onoffEMDm)],[std(offoffEMDm), std(onoffEMDm)],'.k');
xlim([0.5, 2.5]); ylim([0,25]);
ylabel('Wasserstein Distance', 'FontSize', 10)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 10, 'XTick', [1, 2], 'XTickLabels', [{'Control'}, {'Light On'}], 'YTick', 0:10:25);

set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5]);

%Hypothesis Testing
[t, p] = ttest(onoffEMDm,offoffEMDm) % p = 0.0338

%% Figure 4 Heatmaps
%To execute this, load the Summary Electrolytic Lesion Data (produced by electroSumStats.m) file that contains the output data structure (electroStats). The example in the figure uses
%Grn121, which is entry #2 in the input structure.
%This cell code will recalculate the  heatmaps for the conditions, up-sample them, trim to size, plot and format.

%This selects which bird to process
i = 2; 

%Extract the relevant data from the structure
pre = squeeze(electroStats(i).mornN(1,:,:));
postP = squeeze(electroStats(i).postN(5,:,:));
postN = squeeze(electroStats(i).nightN(5,:,:));
post1D = squeeze(electroStats(i).mornN(6,:,:));
post2D = squeeze(electroStats(i).mornN(7,:,:));
post4D =squeeze(electroStats(i).mornN(10,:,:));

%Rescale all for hi-res images
% preBig = imresize(pre, 10);
% postPBig = imresize(postP, 10);
% postNBig = imresize(postN, 10);
% post1DBig = imresize(post1D, 10);
% post2DBig = imresize(post2D, 10);
% post4DBig = imresize(post4D, 10);
preBig = pre;
postPBig = postP;
postNBig = postN;
post1DBig = post1D;
post2DBig = post2D;
post4DBig = post4D;

%Trim to 300ms from 350ms
[r, c] = size(preBig);
subRows = floor(r*(300/350)); %reduce timescale to 300/350
ratio = (subRows/300);

preBig_tr = preBig(1:subRows,:);
postPBig_tr = postPBig(1:subRows,:);
postNBig_tr = postNBig(1:subRows,:);
post1DBig_tr = post1DBig(1:subRows,:);
post2DBig_tr = post2DBig(1:subRows,:);
post4DBig_tr = post4DBig(1:subRows,:);

%Plot to figure
load('inactCmap.mat');
figure(11); clf
subplot(6,1,1); cla
imagesc(preBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
title('Grn121 Heatmaps')

subplot(6,1,2); cla
imagesc(postPBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(6,1,3); cla
imagesc(postNBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(6,1,4); cla
imagesc(post1DBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(6,1,5); cla
imagesc(post2DBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(6,1,6); cla
imagesc(post4DBig_tr', [0, 0.0028]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([1,100,200,300]), 'XTickLabels', [0,100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
xlabel('Syl Dur (ms)'); ylabel('log(Entropy)')

set(gcf, 'Units', 'Inches', 'Position', [0,0,3,10])

%% Figure 3 ElectroLesion Distances Plot
%To execute this, load the Summary Electrolytic Lesion Data (produced by electroSumStats.m) file that contains the output data
%structure (electroStats).
%This cell code will plot the distances between heatmaps in all conditions and format it to be similar to Julie's plots

%Strip out useful info
EMD = [];
for i = 1:length(electroStats)
    EMD(i,:) = electroStats(i).D2emd;
end

%Plot to figure
figure(12); clf
q = [];
b = bar(1:7, mean(EMD,1), 0.3, 'FaceColor', 'none'); hold on
for i = 1:size(EMD,1)
    q(i) = plot(1:7, EMD(i,:), '.k');
    set(q(i), 'Marker', 'o', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'LineStyle', 'none');
end
eb = errorbar(1:7, mean(EMD,1), std(EMD,1)./sqrt(size(EMD,1)),'.k');
xlim([0.5,6.5]); ylim([0,20]); %axis square
ylabel('Wasserstein Distance', 'FontSize', 12)
title(['Summary Inactivation'], 'FontSize', 12);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1, 2, 3, 4, 5, 6, 7], 'XTickLabel', [{'Control'}, {'1st Hour'}, {'8 Hours'}, {'1 Day'}, {'2 Days'}, {'3 Days'} , {'4 Days'}], 'YTick', [0, 20], 'LineWidth', 2, 'FontSize', 10)   

set(gcf, 'Units', 'Inches', 'Position', [0 0 4 2.5])

%Hypothesis Testing
[t.cntrl1st, p.cntrl1st] = ttest(EMD(:,1), EMD(:,2));      % p = 0.0025
[t.cntrl8h, p.cntrl8h] = ttest(EMD(:,1), EMD(:,3));        % p = 0.0162
[t.cntrl1d, p.cntrl1d] = ttest(EMD(:,1), EMD(:,4));        % p = 7.10e-4
[t.cntrl2d, p.cntrl2d] = ttest(EMD(:,1), EMD(:,5));        % p = 0.0558
[t.cntrl3d, p.cntrl3d] = ttest(EMD(:,1), EMD(:,6))         % p = 0.2843

%Do repeated measures one-way ANOVA and run multiple comparisons tests (added 8/11/15)
%Establish repeated measures model
t = table(EMD(:,1), EMD(:,2), EMD(:,3), EMD(:,4), EMD(:,5), EMD(:,6),'VariableNames',{'meas1','meas2','meas3','meas4','meas5','meas6'});
Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas6~1','WithinDesign',Meas);

%Do the RM ANOVA
ranovatbl = ranova(rm);

% Test sphericity assumption using mauchly test -- if p<0.05, sphericity is true and you can use uncorrected pvalue;
% otherwise use corrected values
mauchly(rm)
epsilon(rm)
%                                               SumSq     DF    MeanSq       F         pValue        pValueGG      pValueHF     pValueLB
%                                               ______    __    _______    ______    ________    _________    __________    ________
% 
%     (Intercept):Measurements    38.607     5         7.7214    17.664    8.2028e-06    0.0086922    0.00085415    0.024588
%     Error(Measurements)           6.5568    15    0.43712           
% 
%           W    ChiStat    DF    pValue
%           _    _______    __    ______
%           0          Inf        14    0     
% 
%     Uncorrected    GreenhouseGeisser       HuynhFeldt    LowerBound
%     ___________    _________________    __________    __________
% 
%           1                         0.3303                  0.68311       0.2       



%Do Dunnett's Test with Pre as the control
f = table2array(ranovatbl);
stats = [];
stats.means = [mean(EMD(:,1:6),1)];
stats.df = f(2,2);
stats.n = [size(EMD,1), size(EMD,1), size(EMD,1), size(EMD,1), size(EMD,1), size(EMD,1)];
stats.s = sqrt(f(2,3));
pvals = dunnett(stats, [2:6],1)

% pvals = 
%     Measurements_1    Measurements_2       pValue      
%     ______________    ______________    __________  
% 
%     Cntr                                 Post                3.8771e-06 
%     Cntr                                 8h                   0.0004
%     Cntr                                 1d                   0.0019
%     Cntr                                 2d                   0.0543
%     Cntr                                 3d                   0.1390

    
%% Figure 3 Spontaneous Spiking plot
%To execute this, load the Spontaneous Spiking Data file
%This cell code will plot recovery of spontaneous spiking following lesion

%Get necessary data
% spikes1 = spikeRate2;
% spikes2 = spikeRate3;
% spikes3 = spikeRate4;
% xs = timeZero;
% 
% sz = size(spikes1); s = 0.5;
% spikes1 = spikes1 + (s.*randn(sz)); spikes1(spikes1 < 0) = 0;
% spikes2 = spikes2 + (s.*randn(sz)); spikes2(spikes2 < 0) = 0;
% spikes3 = spikes3 + (s.*randn(sz)); spikes3(spikes3 < 0) = 0;

%Gen Index for day time
% preIdx = timeZero < 0;
% day1Idx = (timeZero >= 0) & (timeZero <= 10.5); 
% day2Idx = (timeZero >= 21.5) & (timeZero <= 34.5);
% 
% %Normalize spike rates by the mean of pre-lesion rate
% spNorm1 = spikes1/(mean(spikes1(preIdx)));
% spNorm2 = spikes2/(mean(spikes2(preIdx)));
% spNorm3 = spikes3/(mean(spikes3(preIdx)));
% 
% %Smooth individual traces with a sliding boxcar; don't over run experiment boundaries
% sp1A = mean(windowTS(spNorm1(day1Idx),5,1,'pad','boxcar')'); 
% sp2A = mean(windowTS(spNorm2(day1Idx),5,1,'pad','boxcar')'); 
% sp3A = mean(windowTS(spNorm3(day1Idx),5,1,'pad','boxcar')'); 
% 
% sp1B = mean(windowTS(spNorm1(day2Idx),7,1,'pad','boxcar')'); 
% sp2B = fliplr(mean(windowTS(spNorm2(day2Idx),7,1,'pad','boxcar')')*.3 + 0.4); 
% sp3B = mean(windowTS(spNorm3(day2Idx),7,1,'pad','boxcar')'); 
% 
% pre(1) = mean(spNorm1(preIdx));
% pre(2) = mean(spNorm2(preIdx));
% pre(3) = mean(spNorm3(preIdx));
% preM = mean(pre); preSEM = std([spNorm1(preIdx), spNorm2(preIdx), spNorm3(preIdx)])/sqrt(3);
% postM = mean([sp1A; sp2A; sp3A]);
% d2M = mean([sp1B; sp2B; sp3B]);
% 
% %Plot
% figure(15); clf
% plot(xs(day1Idx),sp1A, 'Color', [0.5, 0.5, 0.5]); hold on
% plot(xs(day1Idx),sp2A, 'Color', [0.5, 0.5, 0.5])
% plot(xs(day1Idx),sp3A, 'Color', [0.5, 0.5, 0.5])
% plot(xs(day1Idx),postM, 'Color', 'k', 'LineWidth', 3)
% 
% plot(xs(day2Idx),sp1B, 'Color', [0.5, 0.5, 0.5])
% plot(xs(day2Idx),sp2B, 'Color', [0.5, 0.5, 0.5])
% plot(xs(day2Idx),sp3B, 'Color', [0.5, 0.5, 0.5])
% plot(xs(day2Idx),d2M, 'Color', 'k', 'LineWidth', 3)
% b = patch([-2, 34, 34, -2], [1-preSEM, 1-preSEM, 1+preSEM, 1+preSEM], 'k');
% set(b, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.5, 0.5])
% xlim([-1, 33]); ylim([0, 1.3])
% xlabel('Time Post-Lesion'); ylabel('Normalized Firing Rate')
% set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:10:30, 'YTick', 0:0.5:1)
% legend('show'); legend('boxoff')

%% Figure 5 Cengiz's Modeling figures
%These scripts format the data cengiz sent... no idea how his code works or what's required for it, but this hadles the
%formatting for the paper. To execute, load the 'experiment.mat' file he emailed (now in the dropbox folder) and run the cell
%of code

%Load data from file
load experiment_multiple_runs.mat

%Generate Motif Completion Rate Figure
for i = 1:10
    n_chain(i,:) = zeros(1,6900);
    n_chain(i,1:1000) = sum(lengthV(i,1:1000)==80)/1000;

    for j = 1001:6900
        n_chain(i,j) = length(find(lengthV(i,j:min(j+100,7000))==80))/100;
    end
end

chainAm = mean(n_chain(:,1:1000)); chainAs = std(n_chain(:,1:1000))./sqrt(10); 
chainBm = mean(n_chain(:,1001:end)); chainBs = std(n_chain(:,1001:end))./sqrt(10); 

chainAm = mean(windowTS(chainAm, 51, 1, 'pad', 'boxcar')'); chainAs = mean(windowTS(chainAs, 51, 1, 'pad', 'boxcar')');
chainBm = mean(windowTS(chainBm, 51, 1, 'pad', 'boxcar')'); chainBs = mean(windowTS(chainBs, 51, 1, 'pad', 'boxcar')');

figure(1); clf
shadedErrorBar(1:1000, chainAm, chainAs, 'k'); hold on
shadedErrorBar(1001:6900, chainBm, chainBs, 'r')
xlim([0,6900]); ylim([0,1.05]);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',0:2000:6900, 'YTick', [0:0.5:1], 'LineWidth', 1.5, 'FontSize', 16);
xlabel('Number of Renditions', 'FontSize', 16);
ylabel('Normalized Motif Completion Rate', 'FontSize', 16);

set(gcf, 'Units', 'Inches', 'Position',[0,0,5,3])


%Generate Motif Duration plot
for i = 1:10
    d_chain(i,:) = zeros(1,6900);
    index = find(lengthV(i,1:1000)==80);
    d_chain(i,1:1000) = mean(durationV(i,index));

    for j = 1001:6900
        index = find(lengthV(i,j:min(j+100,7000))==80);
        index = index+j-1;
        d_chain(i,j) = mean(durationV(i,index));
    end   
end

normVal = mean(mean(d_chain(:,1:1000)));
chainAm = mean(d_chain(:,1:1000))./normVal; chainAs = std(d_chain(:,1:1000)./normVal)./sqrt(10); 
chainBm = mean(d_chain(:,1001:end))./normVal; chainBs = std(d_chain(:,1001:end)./normVal)./sqrt(10); 

chainAm = mean(windowTS(chainAm, 51, 1, 'pad', 'boxcar')'); chainAs = mean(windowTS(chainAs, 51, 1, 'pad', 'boxcar')');
chainBm = mean(windowTS(chainBm, 51, 1, 'pad', 'boxcar')'); chainBs = mean(windowTS(chainBs, 51, 1, 'pad', 'boxcar')');

figure(2); clf
shadedErrorBar(1:1000, chainAm, chainAs, 'k'); hold on
shadedErrorBar(1001:6900, chainBm, chainBs, 'r')
xlim([0,6900]); ylim([0.99,1.12]);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',0:2000:6900, 'YTick', [1, 1.1], 'LineWidth', 1.5, 'FontSize', 16);
xlabel('Number of Renditions', 'FontSize', 16);
ylabel('Normalized Motif Duration', 'FontSize', 16);

set(gcf, 'Units', 'Inches', 'Position',[0,0,5,3])

%%
%Generate sample data figure
%index_pre = sort(randperm(1000,50));
%index_post = sort(1500+randperm(500,50));
% index_day1 = sort(4000+randperm(1000,50));
%index_day2 = sort(5000+randperm(2000,50));

figure(201); clf
subplot(3,1,1)
durPre = 680*durationV(1,index_pre);
preMean = mean(durPre);
threshFrac = 0.9;
for k = 1:length(durPre)
    b(k) = barh(k, durPre(k), 1, 'EdgeColor', 'none'); hold on
    if durPre(k) > (preMean*threshFrac)
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,425])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:200:1000);

subplot(3,1,2)
durPost = 680*durationV(1,index_post);
for k = 1:length(durPost)
    b(k) = barh(k, durPost(k), 1, 'EdgeColor', 'none'); hold on
    if durPost(k) > (preMean*threshFrac)
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,425])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:200:1000);

subplot(3,1,3)
durDay1 = 680*durationV(1,index_day1);
for k = 1:length(durDay1)
    b(k) = barh(k, durDay1(k), 1, 'EdgeColor', 'none'); hold on
    if durDay1(k) > (preMean*threshFrac)
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,425])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:200:1000);

% subplot(4,1,4)
% durDay2 = 680*durationV(1,index_day2);
% for k = 1:length(durDay2)
%     b(k) = barh(k, durDay2(k), 1, 'EdgeColor', 'none'); hold on
%     if durDay2(k) > (preMean*threshFrac)
%         b(k).FaceColor = [0.5, 0.5, 0.5];
%     else
%         b(k).FaceColor = 'r';
%     end
% end
% axis tight; axis ij
% xlim([0,425])
% set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:200:1000);

%% Figure S2 9nL vs 27nL Inactivation heatmaps 
%To execute this, cd to the home folder (C:\Users\Tim\Desktop\Nif Project Figures\Inactivations\Grn117\)
%and execute section. Script will load the necessary distributions from the file (Pre, and inactivation).
%This cell code will extract the heatmaps for the requested days, trim to size, plot and format.

%Load from file
t = load('Grn117_140820_inactDataset.mat', 'preN', 'injN');
v = load('Grn117_140827_inactDataset.mat', 'injN');
close all

%Copy out data
pre = t.preN;
inj27 = t.injN;
inj9 = v.injN;
clear('t', 'v')
%Filter as required
msize = 7;
sigma = 3;
h = fspecial('gaussian', msize, sigma);

preBig = imfilter(pre,h,'replicate');
inj27Big = imfilter(inj27,h,'replicate');
inj9Big = imfilter(inj9,h,'replicate');

%Limit the range to 0-300ms from the original 0-350ms
[r,c] = size(preBig);
subRows = floor(r*(300/350)); %reduce timescale to 300/350

preBig_tr = preBig(1:subRows,:);
inj27Big_tr = inj27Big(1:subRows,:);
inj9Big_tr = inj9Big(1:subRows,:);

ratio = (subRows/300);

%Plot to figure
load('inactCmap.mat');
figure(12); clf
subplot(3,1,1); cla
imagesc(preBig_tr', [0, 0.0018]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([100,200,300]), 'XTickLabels', [100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])
title('Grn011 Pre-Post Heatmaps')

subplot(3,1,2); cla
imagesc(inj9Big_tr', [0, 0.0018]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([1,100,200,300]), 'XTickLabels', [0,100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

subplot(3,1,3); cla
imagesc(inj27Big_tr', [0, 0.0018]); axis xy; colormap(inactCmap)
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ratio*([1,100,200,300]), 'XTickLabels', [0,100,200,300], 'YTick', [1, c], 'YTickLabels', [-4, 0])

set(gcf, 'Units', 'Inches', 'Position', [0,0,5,4])

%% Plot examples of the HVC recovery and HVC drift for an example bird

%Load the required data from file
folder= 'C:\Users\Tim\Desktop\Temp Nif Data\Grn121\';
lesionData = 'Grn121 HVC Activity (Ch3).mat';
controlData = 'Control121 HVC Activity (Ch3).mat';

L = load([folder lesionData]);
C = load([folder controlData]);

%Samples to average over
numSamp = 25;

%Retrieve mean traces form the dataset
pre = mean(L.data.neuropowAbs((L.data.dayEnds(1)-numSamp):L.data.dayEnds(1),:),1);
post = mean(L.data.neuropowAbs((L.data.dayEnds(1)+1):(L.data.dayEnds(1)+numSamp+1),:),1);
recov = mean(L.data.neuropowAbs(end-numSamp:end,:),1);

drift1 = mean(C.data.neuropowAbs(1:numSamp,:),1);
% drift2 = mean(C.data.neuropowAbs((end-numSamp:end),:),1);

%X-ticks
xpos = 0:4415:length(pre);
xt = 0:100:(length(pre)/44.15);

%Plot figure
figure(13); clf
subplot(1,2,1)
plot(pre); hold on
plot(post);
plot(recov);
axis tight; ylim([0,0.2]);
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Mean HVC Power (V^2)', 'FontSize', 16);
% legend show
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',xpos,'XTickLabel', xt, 'YTick', 0:0.1:0.2, 'LineWidth', 1.5, 'FontSize', 16);

subplot(1,2,2)
plot(drift1); hold on
plot(pre);
axis tight; ylim([0,0.2]);
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Mean HVC Power (V^2)', 'FontSize', 16);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',xpos,'XTickLabel', xt, 'YTick', 0:0.1:0.2, 'LineWidth', 1.5, 'FontSize', 16);

set(gcf, 'Units', 'Inches', 'Position',[0,0,10,4.5])

figure(14); clf
plot(pre, 'DisplayName', 'Pre'); hold on
plot(post, 'DisplayName', 'Post');
plot(recov, 'DisplayName', 'Recov');
plot(drift1, 'DisplayName', 'Drift');
axis tight; ylim([0,0.2]);
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Mean HVC Power (V^2)', 'FontSize', 16);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',xpos,'XTickLabel', xt, 'YTick', 0:0.1:0.2, 'LineWidth', 1.5, 'FontSize', 16);

legend show

set(gcf, 'Units', 'Inches', 'Position',[0,0,5,3])






