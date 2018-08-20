%% Here I plot the full chain propagation rate. Upto run 1000, I plot the 
% overall success rate. After run 1000, I do a moving average of window 
% size 100 and average over 10 simulations.

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


%%
%% Here I plot the full chain duration. Upto run 1000, I plot the 
% average duration. After run 1000, I do a moving average of window 
% size 100, look at only runs for which the chain fully propagated within 
% this window and take a mean among these runs.
% Final plot is the average over 10 simulations.

% clear all;
% close all;
% 
% load experiment_multiple_runs.mat

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
