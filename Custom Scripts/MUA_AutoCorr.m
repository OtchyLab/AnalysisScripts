%Scratch Pad Script to explore autocorrelation lags in MUA recordings
clear

%Specify the files to analyze
motherDir = 'V:\SongbirdData\LW41\Datasets';
fileBase = 'LW41_20180808_dataset_234567ch';

% motherDir = 'C:\Users\Tim\Desktop\Oscillation Analysis';
% fileBase = 'LR81RY177_170830_dataset_ch';

suffix = '.mat';
sets = 1:4;

%Load and preprocess data
xs = [];
for i = sets
    %Select and load file
    file = [motherDir filesep fileBase num2str(i) suffix];
    load(file, 'neuro')
    
    %Preprocess
    mPow = mean(neuro.aligned_abs,1); %Mean over reditions
    mPow_trim = mPow(10:end-10); %Remove the ends (edge effects)
    xs(i,:) = mPow_trim - mean(mPow_trim); %Set mean == 0
    
    %Remove loaded data
    clear('neuro')
end

%Run autocovariance
figure
minLag = 10/1000;
r = []; lags = []; mar = [];
for i = sets
    [r(i,:), l] = xcov(xs(i,:), 'coeff');
    lags = l./44150;
    
    %ID max corr lag
    lim_r = r(lags > minLag);
    lim_lag = lags(lags > minLag);
    [mar(i), idx] = max(abs(lim_r));
    mar_time(i) = lim_lag(idx);
    
    %Plot output
    plot(lags, r(i,:)); hold on
end
m_lags = lags;
m_r = mean(r,1);
plot(m_lags, m_r, 'r', 'LineWidth', 2)
xlim([0,.5])
xlabel lags; ylabel correlation;

%Identify the lag of maximum correlation
lim_r = m_r(m_lags > minLag);
lim_lag = m_lags(m_lags > minLag);
[m, idx] = max(abs(lim_r));
m_time = lim_lag(idx);

a = 1;



% for i = 1 : 100000
% xb = x(randperm(numel(x)));
% [rb, lags] = xcov(xb, 'coeff');
% marb(i) = max(abs(rb(lags > 0)));
% end
% pval = mean(marb >= mar)
% 
% rc = quantile(marb, 0.95)

%%

file = dir('*1.mat');
traces = [];
figure(10000); clf
for i = 1:numel(file)
    load(file(i).name, 'neuro')
    
    traces(i,:) = 1000.*mean(neuro.aligned_abs(:,10:end-10),1);
%     traces(i,:) = 1000.*mean(neuro.aligned(:,10:end-10),1);

    plot((1:numel(traces(i,:)))./44150, traces(i,:), 'LineWidth', 1.5); hold on
end
axis tight; %ylim([0, 0.5])
xlabel('Time (ms)'); ylabel('V^2 (mv^2)')
set(gca, 'Box', 'off', 'TickDir', 'out')


for i = 1:numel(file)
    rho(i) = corr(traces(1,:)', traces(i,:)');

    mAmp(i) = mean(traces(i,:));

end
%%
figure(10001); clf
subplot(1,2,1)
xs = [1, 5, 10, 15];
plot(xs, rho, '-ok', 'LineWidth', 2)
xlim([0,16]); ylim([0,1])
xlabel('Time (DPI)'); ylabel('R^2')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(1,2,2)
plot(xs, mAmp, '-ok', 'LineWidth', 2)
xlim([0,16]); ylim([0,0.2])
xlabel('Time (DPI)'); ylabel('Mean V^2')
set(gca, 'Box', 'off', 'TickDir', 'out')





