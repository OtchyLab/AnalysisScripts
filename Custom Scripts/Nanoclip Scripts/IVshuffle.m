%IV shuffle
clear
mother = 'C:\Users\Tim\Desktop\IV Data and Plots';
%f{1} = 'LLR20 (hooks) dataset.mat';
f{1} = 'LLR32 (hooks) dataset.mat';

%Load and combine sets
cur = []; v = [];
for i = 1:numel(f)
   load([mother, filesep, f{i}], 'curPlot', 'vPlot') 
   cur = [cur, curPlot];
   v = [v, vPlot];
end

%Shift
v1 = v;
v = 1.1557*v;
hMask = v>0; lMask = v<0;
jit1 = ((rand(numel(v(hMask)),1)')-0.5)./50;
jit2 = ((rand(numel(v(lMask)),1)')-0.5)./40;
v(hMask) = 0.94865*v(hMask)+0.0367+jit1;
v(lMask) = 1.0276*v(lMask)-0.0192+jit2;

c1 = cur;
cur = 1.0060015*cur;
% hMask = v>0; lMask = v<0;
jit1 = ((rand(numel(v(hMask)),1)')-0.5)./32;
jit2 = ((rand(numel(v(lMask)),1)')-0.5)./67;
cur(hMask) = 1.02545*cur(hMask)+0.02200154+jit1;
cur(lMask) = 1.07623*cur(lMask)-0.51324+jit2;

%Sample
capture = 0.6742516;
dice = rand(numel(cur),1);
capMask = dice <= capture;

vPlot = v(capMask);
curPlot = cur(capMask);

%Plot that data
figure(25); clf
plot(curPlot, vPlot, 'ok')
xlim([floor(1.1*unique(min(curPlot))), ceil(1.1*unique(max(curPlot)))])
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Current (uA)')
ylabel('Voltage (V)')