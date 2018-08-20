function mockDistShifts(startMean,startStd,endMean,endStd)

%Create pool of starting data
n = 100;
init_rand = (randn(n,1) + (rand(n,1))); %Normal distribution + white noise
startPool = (init_rand*startStd)+startMean;
[startCounts, startBins] = hist(startPool,540:5:620);

m = 136;
end_rand = (randn(m,1) + (rand(m,1))); %Normal distribution + white noise
endPool = (end_rand*endStd)+endMean;
[endCounts, endBins] = hist(endPool,540:5:620);

%Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');

%Blue distribution plot
bar(startBins(startCounts~=0),startCounts(startCounts~=0),'FaceColor','none','EdgeColor',[0 0 1]);

%Red distribution plot
bar(endBins(endCounts~=0),endCounts(endCounts~=0),'FaceColor','none','EdgeColor',[1 0 0]);

%Plot dashed lines for the means of each distribution
maxy = get(gca,'YLim');
sMean = mean(startPool);
eMean = mean(endPool);
line([sMean,sMean],[0,maxy(2)],'Color',[0 0 0],'LineStyle','--');
line([eMean,eMean],[0,maxy(2)],'Color',[0,0,0],'LineStyle','--');
set(gca,'TickDir','out')