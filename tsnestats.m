%Process figures and descriptive stats for the nanocliip fictive singing TSNE plots

%Reset workspace
clear

%Load the raw data
load('/Users/Tim/Desktop/tsne_summary_data')

%Initialize output structure
stats = [];

numBirds = numel(tsne_summary);
within = getFieldVector(tsne_summary,'within');
across = getFieldVector(tsne_summary,'across');
data = [within', across'];

%Calculate descriptive stats
%Mean
stats.mean = mean(data,1);

%Median
stats.median = median(data,1);

%Std
stats.std = std(data,1);

%SEM
stats.sem = stats.std ./ sqrt(numBirds);

%Significance testing
[h_norm_within, p_norm_within] = kstest(data(:,1)); %[h_norm_within, p_norm_within] = [1, 0]
[h_norm_across, p_norm_across] = kstest(data(:,2)); %[h_norm_across, p_norm_across] = [1, 8e-16]

[h_signif, p_signif] = ttest(data(:,1), data(:,2)); %[h_signif, p_signif] = [1, 0.0022]

%Plot output
figure(69); clf
set(gcf,'Units', 'Inches', 'Position', [15, 4, 4, 5.75]);

b = bar([1,2], stats.mean); hold on
eb = errorbar([1,2], stats.mean, stats.sem, '.k');

xlim([0.5,2.5])
ylim([0, 40])
ylabel('Mean Distance (AU)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [1,2], 'XTickLabel', {'Within'; 'Across'})
set(gca, 'YTick', [0, 40], 'FontSize', 12, 'LineWidth', 1.5)




