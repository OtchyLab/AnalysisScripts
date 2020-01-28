%IV shuffle4
clear
mother = 'C:\Users\Tim\Desktop\IV Data and Plots';
f = 'All Base (nc) dataset.mat';

%Load and combine sets
cur = []; v = [];
load([mother, filesep, f], 'curPlot', 'vPlot')
cur = curPlot;
v = vPlot;

emaskL = (cur > -34.9 & cur < -34.5) | (cur > -29.84 & cur < -29.32) | (cur > -24.95 & cur < -24.50) | (cur > -20.03 & cur < -19.7) | (cur > -16.7 & cur < -16.45) | (cur > -16.09 & cur < -15.85) | (cur > -15 & cur < -14.85) | (cur > -9.95 & cur < -9.74) | (cur > -7.9 & cur < -7.8)| (cur > -7.4 & cur < -7.3);
emaskH = (cur > 4.9 & cur < 5) | (cur > 5.9 & cur < 6) | (cur > 7.45 & cur < 7.6) | (cur > 15.95 & cur < 16.2) | (cur > 16.65 & cur < 16.8) | (cur > 19.85 & cur < 20.03) | (cur > 24.8 & cur < 25.02) | (cur > 34.6 & cur < 35) | (cur > 46.7 & cur < 47);
emaskTot = emaskL | emaskH;
cur = cur(~emaskTot);
v = v(~emaskTot);

%Shift
v1 = v;
v = 0.52*v;
hMask = v>0; lMask = v<0;
jit1 = ((rand(numel(v(hMask)),1)')-0.5)./10;
jit2 = ((rand(numel(v(lMask)),1)')-0.5)./10;
v(hMask) = 1.0578*v(hMask)-0.207+jit1;
v(lMask) = 0.87545*v(lMask)+0.0233+jit2;

c1 = cur;
cur = 0.66556*cur;
% hMask = v>0; lMask = v<0;
jit1 = ((rand(numel(v(hMask)),1)')-0.5)./15;
jit2 = ((rand(numel(v(lMask)),1)')-0.5)./7;
cur(hMask) = 1.04545*cur(hMask)+0.07154+jit1;
cur(lMask) = 1.01623*cur(lMask)+0.071324+jit2;

%Sample
capture = 0.22;
dice = rand(numel(cur),1);
capMask = dice <= capture;

vPlot = v(capMask);
curPlot = cur(capMask);

%Plot that data
figure(31); clf
plot(curPlot(37:end-21), vPlot(37:end-21), 'ok')
xlim([floor(1.1*unique(min(curPlot))), ceil(1.1*unique(max(curPlot)))])
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Current (uA)')
ylabel('Voltage (V)')