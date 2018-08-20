function xcorrSummaryOverBirdsPlots()

%Get the Processed file to analyze
[temp,path] = uigetfile('*.mat','Select the processed data file.','MultiSelect','on');

if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames.name = temp;
end

summary_lags = [];
summary_means = [];
summary_zscores = [];
summary_z_stack = [];
new_zscores = [];

fig1 = figure;
hold on
%Load everything from file, process, and plot
for i = 1:length(fnames)
    %Load the data
    load([path filesep fnames(i).name],'lag_neuroCov','m_zstack')
    summary_lags = [summary_lags;lag_neuroCov];
    summary_zscores = [summary_zscores;m_zstack];
    
    %Covert the mean into a z-score
    z = (m_zstack-mean(m_zstack))./std(m_zstack);
    new_zscores = [new_zscores;z];
    
    %Show data on axes
    h = scatter(lag_neuroCov,z,'Marker','.');
end
title(['Scatter Plot of all Days/Birds'])
xlabel('Neural Signal Lead (ms)')
ylabel('Z-Score (x-m)/s')
set(gca,'TickDir','out','Box','off')

%Create summary average plot
summary_m_zstack = mean(new_zscores,1);
summary_sem_zstack = std(new_zscores,1)/sqrt(i);

fig3 = figure;
hold on
%scatter(lags(1,:),mean(zscores,1),'Marker','x')
errorbar(lag_neuroCov,summary_m_zstack,summary_sem_zstack,'s')
title(['Lags overall birds and days'])
xlabel('Neural Signal Lead (ms)')
ylabel('Z-Score (x-m)/s')
set(gca,'TickDir','out','Box','off')

%xs = -36:.05:100;
xs = -36:.05:100;
interped = interp1(lag_neuroCov,summary_m_zstack,xs);
% sm_series = smooth(interped,5,'moving');
% plot(xs,sm_series,'r')
% [~,I] = max(sm_series);
% maxLagSMOOTH = xs(I)

[sigma,mu,A]=mygaussfit(xs,(interped+abs(min(interped))));
y=A*exp(-(xs-mu).^2/(2*sigma^2));
plot(xs,y-abs(min(interped)),'g')
[~,i] = max(y);
maxLagGauss = xs(i)
 
%Save image and data to file
save([path 'Summary Data Over Birds and Days.mat'],'lag_neuroCov','summary_m_zstack','summary_sem_zstack','fnames','maxLagGauss')
saveas(fig1,[path 'Summary Data Over Birds and Days Breakout.fig'])
saveas(fig3,[path 'Summary Data Over Birds and Days.fig'])


