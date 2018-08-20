
figure(100); clf

%Plot clip birds
block = 1:3;
b = [];
for i = block
%     plot(summaryData(i).days, summaryData(i).MeanSpect, ':.b', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
mS = mean(mean(b(1:3,:), 2));
shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'b', 1); hold on

%Plot isolation birds
block = 4:6;
b = [];
for i = block
%     plot(summaryData(i).days, summaryData(i).MeanSpect, ':.k', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
mS = mean(mean(b(1:3,:), 2));
shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'k', 1)

%Plot crush birds
block = 7:9;
b = [];
for i = block
%     plot(summaryData(i).days, summaryData(i).MeanSpect, ':.r', 'LineWidth', 1); hold on
    b = [b, summaryData(i).MeanSpect];
end
mS = mean(mean(b(1:3,:), 2));
shadedErrorBar(summaryData(i).days, mean(b, 2)./mS, std(b, 1, 2)./mS, 'r', 1)

%Format
xlim([-3.5, 8.5]); ylim([0.7, 1.05])
xlabel('Time (days)'); ylabel('Acoustic Similarity')
set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0.6:0.1:1)

set(gcf, 'Units', 'Inches', 'Position', [1, 1, 4, 4])





