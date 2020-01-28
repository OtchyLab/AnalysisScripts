%IVtrim
% %Plot data
% figure(29); clf
% [edgesX2,edgesY2,N,histImgHandle] = ndhist(curPlot, vPlot, 'bins', 10, 'filter', 'max');

%Draw ROI to remove
ROI = imfreehand(gca);
ROIposition = getPosition(ROI);

%Translate the position of the ROI into image coordinates
%Translate the position of the ROI into image coordinates
for i = 1:size(ROIposition,1)
    [~, indx(i)] = min(abs((edgesX2-ROIposition(i,1))));
    [~, indy(i)] = min(abs((edgesY2-ROIposition(i,2))));
end

%Mask gen
mask = poly2mask(indx, indy, length(edgesY2), length(edgesX2));

%Generate index of points that fall within the mask
selInd = [];
for i =1:length(curPlot)
    [~, pntsx(i)] = min(abs((edgesX2-curPlot(i))));
    [~, pntsy(i)] = min(abs((edgesY2-vPlot(i))));
    selInd(i) = mask(pntsy(i),pntsx(i));
end

% newCur = curPlot(~selTot');
% newV = vPlot(~selTot');
% figure(30); clf
% [edgesX2,edgesY2,N,histImgHandle] = ndhist(newCur, newV, 'bins', 10);
