%Make the spontaneous activity example figure


%Filter data 
%Set filtering params
fs = 44150;
HP_fNorm = 1000/(fs/2); % 1kHz
LP_fNorm = 5000/(fs/2); % 15kHz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]); %Butterworth, 2-p

preFilt = filtfilt(BP_b,BP_a, pre);
postImFilt = filtfilt(BP_b,BP_a, postImmediate);
post1Filt = filtfilt(BP_b,BP_a, post1);
post3Filt = filtfilt(BP_b,BP_a, post3);
post24Filt = filtfilt(BP_b,BP_a, post24);



%Set indices for equal-length plotting
preIndx = (1:300001);
postImIndx = (1:300001);
post1Indx = (5530:305530);
post3Indx = (335089:635089);
post24Indx = (246478:546478);

ticks= (0:6);

figure;
subplot(5,1,1)
plot(preFilt(preIndx))
axis tight; ylim([-0.4, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:fs:length(preIndx), 'XTickLabel', [])
ylabel('Pre')

subplot(5,1,2)
plot(postImFilt(postImIndx))
axis tight; ylim([-0.4, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:fs:length(preIndx), 'XTickLabel', [])
ylabel('Immediate')

subplot(5,1,3)
plot(post1Filt(postIndx))
axis tight; ylim([-0.4, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:fs:length(preIndx), 'XTickLabel', [])
ylabel('Post 1h')

subplot(5,1,4)
plot(post3Filt(post3Indx))
axis tight; ylim([-0.4, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:fs:length(preIndx), 'XTickLabel', [])
ylabel('Post 3h')

subplot(5,1,5)
plot(post24Filt(post24Indx))
axis tight; ylim([-0.4, 0.4])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 1:fs:length(preIndx), 'XTickLabel', ticks);
xlabel('Time (s)')
ylabel('Post 1d')

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4, 6])
