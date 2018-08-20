function sponSpikeSumStats(sponSumStats)

%Define outputs
pre = nan(length(sponSumStats),8);
post = nan(length(sponSumStats),41);
post1 = nan(length(sponSumStats),51);

for i = 1:length(sponSumStats)
    %Extract pre
    t = sponSumStats(i).preRate;
    startPos = 8 - length(t)+ 1;
    pre(i,startPos:end) = t;
    
    %Extract post
    t = sponSumStats(i).postRate;
%     startPos = 41 - length(t)+ 1;
    post(i,1:length(t)) = t;
    
    %Recover time constants from smoothed recovery
    r = nanmean(windowTS(t, 5, 1, 'pad', 'boxcar'), 2);
    range = max(r) - min(r);
    id = find(r>(0.63*range + min(r)), 1, 'first')
    tau(i) = id*(15/60); %calculate tau in hours
    
    %Extract post
    t = sponSumStats(i).post1Rate;
%     startPos = 51 - length(t)+ 1;
    post1(i,1:length(t)) = t;
    
end

%Recover time constants
tau_m = mean(tau);
tau_std = std(tau)

%Across bird stats
m_pre = nanmean(pre,1);
m_post = nanmean(post,1);
m_post1 = nanmean(post1,1);
s_pre = nanstd(pre,1);
s_post = nanstd(post,1);
s_post1 = nanstd(post1,1);

%smoothed means
sm_pre = nanmean(windowTS(m_pre, 5, 1, 'pad', 'boxcar'), 2);
sm_post = nanmean(windowTS(m_post, 5, 1, 'pad', 'boxcar'), 2);
sm_post1 = nanmean(windowTS(m_post1, 5, 1, 'pad', 'boxcar'), 2);

%smoothed sem
ss_pre = nanmean(windowTS(s_pre, 5, 1, 'pad', 'boxcar'), 2)/sqrt(length(sponSumStats));
ss_post = nanmean(windowTS(s_post, 5, 1, 'pad', 'boxcar'), 2)/sqrt(length(sponSumStats));
ss_post1 = nanmean(windowTS(s_post1, 5, 1, 'pad', 'boxcar'), 2);

%Plotting times
preTime = -2:0.25:-0.25;
postTime = 0:0.25:10;
post1Time = 12:0.25:24.5;

%Plot output
figure(1); clf
shadedErrorBar(preTime, sm_pre, ss_pre, 'k'); hold on
shadedErrorBar(postTime, sm_post, ss_post, 'r');
shadedErrorBar(post1Time(1:end-2), sm_post1(1:end-2), ss_post1(1:end-2), 'r');
xlim([-2, 24.5]); ylim([0, 1.3])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', sort([0,-2:4:26]), 'YTick', 0:0.5:1)
xlabel('Time Post Lesion (h)')
ylabel('Normalized Spiking Rate')

set(gcf, 'Units', 'Inches', 'Position', [0,0,4, 4])


