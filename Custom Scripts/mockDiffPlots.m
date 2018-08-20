function mockDiffPlots

%Create figure
figure
hold on

%Generate data
n = 4;
pool = randn(1000,1);
%x = pool(find(pool>.72 & pool<.86,n));
x = pool(find(pool>.25 & pool<.6,n));

pool = randn(1000,1);
y = pool(find(pool>.77 & pool<.91,n));

pool = randn(1000,1);
z = pool(find(pool>.86 & pool<.96,n));

% %Plot data
% plot([0,1],[0,1],'k') %unity line
% hold on
% scatter(x,y)
% 
% %Format plots
% axis tight
% set(gca,'XTick',[0,1],'YTick',[0,1])
% set(gca,'TickDir','out','Box','off')

target_unalignedM = mean(x)
target_unalignedS = std(x)

target_alignedM = mean(y)
target_alignedS = std(y)

nontarget_alignedM = mean(z)
nontarget_unalignedS = std(z)

%Plot the mean and std bars
bar([0],mean(x),'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar([1],mean(y),'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar([2],mean(z),'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);

errorbar(0,mean(x),std(x),'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(1,mean(y),std(y),'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(2,mean(z),std(z),'Marker','.','LineStyle','none','Color',[0 0 0])

%Plot the individual points on top
for i = 1:n
    plot([0,1,2],[x(i), y(i), z(i)],'-ok')
end

%Format plots
axis tight
xlim([-.5,2.5])
ylim([0,1])
set(gca,'XTick', [],'YTick',[0,.5,1])
set(gca,'TickDir','out','Box','off')



