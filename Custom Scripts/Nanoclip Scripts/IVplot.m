%IVplot
clear
mother = 'C:\Users\Tim\Desktop\IV Data and Plots';
f{1} = 'LLR20 (hooks) dataset.mat';
f{2} = 'LLR32 (hooks) dataset.mat';
f{3} = 'LLY72 (hooks) dataset.mat';

fn{1} = 'LLBlk07 (nc) dataset.mat';
fn{2} = 'LLBlk53 (nc) dataset.mat';
fn{3} = 'LLR32 (nc) dataset.mat';

%Load and plot sets
cur = []; v = [];
figure(68); clf
for i = 1:numel(f)
   load([mother, filesep, fn{i}], 'curPlot', 'vPlot') 
   cur = [cur, curPlot];
   v = [v, vPlot];
   
   plot(curPlot, vPlot, 'o'); hold on
end
xlim([floor(1.1*unique(min(curPlot))), ceil(1.1*unique(max(curPlot)))])
set(gca, 'Box', 'off', 'TickDir', 'out')
title('IV Curves for Nanoclip Electrodes (n=3)')
xlabel('Current (uA)')
ylabel('Voltage (V)')

%Calculate the mean and std over devices
% [curS, idx] = sort(cur);
% vS = v(idx);
