function PlotReachAttribClust(OrderedFiles, AttributeImage, ClusterSet)

%Initialize figure
figure;

%Plot the Reachability Plot on upper subplot
subplot(2,1,1);
bar(OrderedFiles(:,3), 'b');
xlim([1 length(OrderedFiles)]);
ylim([0 round(max(OrderedFiles(:,3))+1.5)]);
box on;
hold on;

numClust = length(ClusterSet);
spacing = 7:(1/numClust):8;
for i=1:numClust
    lengthbar = length(ClusterSet(i,1):ClusterSet(i,2));
    plotbar = ones(lengthbar,1)*spacing(i);
    plot(ClusterSet(i,1):ClusterSet(i,2), plotbar, '-r');
    hold on
end
% Create title
title('Reachability Plot fo PCASet, Eps = [], MntPnts=200');


%Plot the Attribute Image on lower subplot
subplot(2,1,2)
imagesc(AttributeImage,[0 255]); colormap(gray);
xlim([1 length(OrderedFiles)]);
ylim([.5 size(AttributeImage,1)+.5]);
box on;
hold on;

% Create xlabel
xlabel('Ordered Points (n-->)');

end %temp end function