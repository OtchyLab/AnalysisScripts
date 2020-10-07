function acute_sweep_summary
%This function will take as input the output of function acute_sweep_analysis_thinCMS.m, 
%and create a summary plot and across-animal statistics. Figures to be used in
%Rowan et al (2020) paper.
%
%
% Updated by TMO 10/05/20

%Load the summary stats from file
folder = '/Users/tim/Desktop/Acute Recordings';
file = 'Sweep Summary Data.mat';
load([folder, filesep, file])

%Plot the summary figure
f = figure(1); clf
sym = {'o', 'x', 'sq', 'd'};
set(gcf, 'Units', 'Inches', 'Position', [10, 10, 5.5, 3.5])

for i = 1:size(binned_out, 1)
    jit = (rand(1,numel(centers))-0.5);
    plot(centers + jit, binned_out(i,:), [sym{i}], 'Color', [0.75, 0.75, 0.75], 'MarkerSize', 9); hold on
end

errorbar(centers, mean(binned_out, 1), std(binned_out, 1), '-k')

xlim([0, 110])
ylim([0, .75])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:50:100, 'YTick', [0, 0.75])
xlabel('Current (\muA)')
ylabel('Vpp (mV)')

