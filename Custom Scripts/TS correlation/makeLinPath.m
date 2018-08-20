function [linpath] = makeLinPath(tbN1,tbN2,tb1)

%Point shift between first and second shift
offsets = tbN2-tbN1;

%Cross points are simply the interpolation between the two closest points
%to the syllable marker
%crossPnts = interp1(tbN1,offsets,tb1,'linear','extrap')';

%Cross points are the weighted average of all markers within a window
%around the syllable marker; all weights outside the window are 0
win = 20; %1/2 window size, in ms
for i = 1:length(tb1)
    dist = abs(tbN1-tb1(i));
    weights = ((-1/win)*dist)+1;
    mask = weights<0;
    weights(mask) = 0;
    crossPnts(i) = sum(weights.*offsets)/sum(weights);
end

%Convert back to path function
linpath = crossPnts'+tb1;

%Plot the output
figure
scatter(tbN1,offsets);hold on
ys = ones(length(tb1),1)*[min(offsets),max(offsets)];
line([tb1,tb1]',ys','Color','k','LineStyle','--','LineWidth',1.5)
scatter(tb1,crossPnts,100,'sr')
plot(tb1,crossPnts,'r')
axis tight

