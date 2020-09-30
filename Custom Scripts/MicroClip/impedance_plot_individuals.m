function impedance_plot_individuals
% This function will take as input a previously created source file containing
% the impedances of each electrode in a set of microclips, do some summary 
% stats on the values, and plot the results. Easy peasy.
%
%Created by TMO 09/30/20

%Source data
folder = '/Users/tim/Desktop/Acute Recordings';
file = 'impedances.mat';

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

%Plot simple figure with means for each device
h1 = plot_simple_means(raw, plated, invivo, post);

%Plot figure with individual points and device means
h2 = plot_points_and_means(raw, plated, invivo, post);

function f = plot_simple_means(raw, plated, invivo, post)
%Initialize figure
f = figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [2.4, 5.6, 5, 5.5])

%Plotting colors
num_dev = size(raw,1);
col = distinguishable_colors(num_dev); %indicates device; RGB triplet

%Plot the means over pads for each device and condition
for i = 1:num_dev
    jit = (2 * rand(1) - 1) * 0.08;
    scatter([1:4] + jit, [nanmean(raw(i,:)), nanmean(plated(i,:)), nanmean(invivo(i,:)), nanmean(post(i,:))], 36, col(i,:), 'filled');
    hold on
end

%Format the figure per usual
format_fig(gca)

function f = plot_points_and_means(raw, plated, invivo, post)
%Initialize figure
f = figure(2); clf
set(gcf, 'Units', 'Inches', 'Position', [7.75, 5.6, 5, 5.5])

%Plotting colors
num_dev = size(raw,1);
num_pad = size(raw,2);
col = distinguishable_colors(num_dev); %indicates device; RGB triplet

%Cycle over devices
for i = 1:num_dev
    %Plot the device means
    jit = (2 * rand(1) - 1) * 0.25;
    S = scatter([1:4] + jit, [nanmean(raw(i,:)), nanmean(plated(i,:)), nanmean(invivo(i,:)), nanmean(post(i,:))], 64, [0.75, 0.75, 0.75], 'filled'); %'MarkerEdgeColor', col(i,:), 'filled');
    S.MarkerEdgeColor = col(i,:);
    hold on
    
   %for j = 1:num_pad 
        %Plot individual pad points
        jit_i = (2 * rand(1,num_pad) - 1) * 0.1;
        plot(1*ones(1,num_pad)+jit+jit_i, raw(i,:), '.', 'Color', col(i,:))
        plot(2*ones(1,num_pad)+jit+jit_i, plated(i,:), '.', 'Color', col(i,:))
        plot(3*ones(1,num_pad)+jit+jit_i, invivo(i,:), '.', 'Color', col(i,:))
        plot(4*ones(1,num_pad)+jit+jit_i, post(i,:), '.', 'Color', col(i,:))
   %end
end

%Format the figure per usual
format_fig(gca)

function format_fig(ax)
%Limits
xlim([0.5, 4.5]); ylim([0, 1000])

%Format
set(ax, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:500:1000)
set(ax, 'XTick', 1:4, 'XTickLabels', {'Raw', 'Plated', 'In-Vivo', 'Post'})
ylabel('Impedance (k\Omega)')



