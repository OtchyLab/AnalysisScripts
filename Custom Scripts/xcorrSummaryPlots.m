function xcorrSummaryPlots()

%Get the Processed file to analyze
[temp,path] = uigetfile('*.mat','Select the processed data file.','MultiSelect','on');

if size(temp,2)>1
    for i = 1:size(temp,2)
        fnames(i).name = char(temp{i});
    end
else
    fnames.name = temp;
end

lags = [];
means = [];
stds = [];
names = [];
zscores = [];
z_stack = [];

fig1 = figure;
hold on
%Load everything from file, process, and plot
for i = 1:length(fnames)
    %Load the data
    load([path filesep fnames(i).name])
    names = [names,name_neuroCov];
    lags = [lags,lag_neuroCov];
    stds = [stds,std_neuroCov];
    means = [means,m_neuroCov];
    
    %Covert the mean into a z-score
    z = (m_neuroCov-mean(m_neuroCov))./std(m_neuroCov);
    zscores = [zscores,z];
    z_stack = [z_stack;z];
    
    %Show data on axes
    h = scatter(lag_neuroCov,z,'Marker','.','DisplayName',name_neuroCov{1}(1:17));
end
title(['Scatter Plot of all Days/Birds'])
xlabel('Neural Signal Lead (ms)')
ylabel('Z-Score (x-m)/s')
set(gca,'TickDir','out','Box','off')

% %Create summary plot
% fig2 = figure;
% scatter(lags,zscores,'Marker','.')
% title(['Scatter Plot of all Days/Birds'])
% xlabel('Neural Signal Lag (ms)')
% ylabel('Z-Score (x-m)/s')

%Create summary average plot
m_zstack = mean(z_stack,1);
s_zstack = std(z_stack,1);

fig3 = figure;
hold on
%scatter(lags(1,:),mean(zscores,1),'Marker','x')
errorbar(lag_neuroCov,m_zstack,s_zstack,'x')
title(['Mean Fits Over 4 Days (Pur626)'])
xlabel('Neural Signal Lead (ms)')
ylabel('Z-Score (x-m)/s')
set(gca,'TickDir','out','Box','off')

xs = -36:4:100;
interped = interp1(lag_neuroCov,m_zstack,xs);
sm_series = smooth(interped,5,'moving');
plot(xs,sm_series,'r')
[~,I] = max(sm_series);
maxLagSMOOTH = xs(I);

[sigma,mu,A]=mygaussfit(xs,(interped+abs(min(interped))));
y=A*exp(-(xs-mu).^2/(2*sigma^2));
plot(xs,y-abs(min(interped)),'g')
[~,i] = max(y);
maxLagGauss = xs(i);

% bird = 'Pur626';
% %Save image and data to file
% save([path bird ' Summary Data.mat'],'names','lags','stds','means','zscores','lag_neuroCov','m_zstack')
% saveas(fig1,[path bird ' Summary Data Breakout.fig'])
% saveas(fig3,[path bird ' Summary Data.fig'])


