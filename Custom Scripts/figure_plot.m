% clear all;
% close all;
% 
% load experiment.mat

n_chain = zeros(1,6900);
n_chain(1:1000) = 1;

for i = 1001:6900
   
    n_chain(i) = length(find(chain_length(i:min(i+100,7000))==80))/100;
    
end
n_chain = 2-n_chain;

%Do smoothing on the original data
chainB = mean(windowTS(n_chain(1001:6900), 51, 1, 'pad', 'boxcar')');

%Plot song stretching figure
figure(1); clf
% subplot(2,1,1); cla; %original
% plot(1:6900,n_chain,'-'); 
% axis([0 7000 1 2])
% 
% subplot(2,1,2); cla
plot(1:1000, n_chain(1:1000), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Pre-Lesion'); hold on
plot(1001:6900, chainB, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Post-Lesion');
xlim([0,6900]); ylim([0.98,2]);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',0:2000:6900, 'YTick', [1,1.75], 'LineWidth', 1.5, 'FontSize', 16);
xlabel('Number of Renditions', 'FontSize', 16);
ylabel('Normalized Motif Truncation Rate', 'FontSize', 16);

set(gcf, 'Units', 'Inches', 'Position',[0,0,5,3])

% exportfig(gcf,'./success_rate.eps','height',1,'width',5,'Color','rgb')

ind = find(chain_length == 80);

t = chain_duration(ind)/mean(chain_duration(1:1000));
tA = t(1:1000);
tB = t(1001:end);
durA_sm = mean(windowTS(tA, 51, 1, 'pad', 'boxcar')'); 
durB_sm = mean(windowTS(tB, 51, 1, 'pad', 'boxcar')');

figure(2); clf
% subplot(2,1,1); cla
% plot(ind,chain_duration(ind)/mean(chain_duration(1:1000)),'.') %original
% % exportfig(gcf,'./duration.eps','height',1,'width',5,'Color','rgb')
% 
% subplot(2,1,2); cla
plot(ind(1:1000), durA_sm, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'DisplayName', 'Pre-Lesion'); hold on
plot(ind(1001:(length(durB_sm)+1000)), durB_sm, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Post-Lesion');
xlim([0,6900]); ylim([0.99,1.12]);
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick',0:2000:6900, 'YTick', [1,1.10], 'LineWidth', 1.5, 'FontSize', 16)
xlabel('Number of Renditions', 'FontSize', 16)
ylabel('Normalized Motif Stretch', 'FontSize', 16)

set(gcf, 'Units', 'Inches', 'Position',[0,0,5,3])


%%

% plot(Thdist(:,10))
% axis([0 1201 -0.051 -0.0499])
% exportfig(gcf,'./pre_th.eps','height',1,'width',2,'Color','rgb')
% 
% plot(Thdist(:,end))
% axis([0 1201 -0.051 -0.0499])
% exportfig(gcf,'./post_th.eps','height',1,'width',2,'Color','rgb')

