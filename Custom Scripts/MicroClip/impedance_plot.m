function impedance_plot
%This function will take as input a previously created source file containing
% the impedances of each electrode in a set of wrapClipsa, do some summary 
% stats on the values, and plot the results. Easy peasy.
%
%Created by TMO 08/26/20

%Source data
folder = '/Users/tim/Desktop/Acute Recordings';
file = 'impedances.mat';

%Plotting options
sym = {'o', 's'}; %indicates electrode layer
col = {'b', 'r', 'g', 'k', 'm', 'c', 'y'}; %indicates device


%Load data from file
load([folder, filesep, file])
if ~exist('impedance')
    disp('Oops... missing the expected data. Have a look and retry.')
    return
end

%Extract fields
raw = cell2mat(getFieldVectorCell(impedance, 'raw'))';
plated = cell2mat(getFieldVectorCell(impedance, 'plated'))';
invivo = cell2mat(getFieldVectorCell(impedance, 'invivo'))';
post = cell2mat(getFieldVectorCell(impedance, 'post'))';

%Simple mean plot
figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [2.4    5.6    5    5.5])
for i = 1:numel(impedance)
    jit = (2 * rand(1) - 1) * 0.08;
    scatter([1:4] + jit, [nanmean(raw(i,:)), nanmean(plated(i,:)), nanmean(invivo(i,:)), nanmean(post(i,:))], ['o', col{i}], 'filled');
    hold on
end

xlim([0.5, 4.5]); ylim([0, 600])
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:250:500)
set(gca, 'XTick', 1:4, 'XTickLabels', {'Raw', 'Plated', 'In-Vivo', 'Post'})
ylabel('Impedance (k\Omega)')



