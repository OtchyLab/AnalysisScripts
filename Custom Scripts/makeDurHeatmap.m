function [matrix, upperScale] = makeDurHeatmap(TargetDur)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot heatmap of durations
%           Based on Risa's code
%           TargetDur = a vector whose elements are the units to be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter parameters
window = 7;
factor = 6;

%Other parameters
maxInterval = 400;
trialBinSize = 50;
intertapBinSize = 5;
upperScale = 7.5E-2;

intertapMatrix = zeros(length(TargetDur), maxInterval);
convMatrix = ones(trialBinSize, intertapBinSize);
for m = 1:length(TargetDur)
    if (TargetDur(m) > 0) && (TargetDur(m) < maxInterval)
        intertapMatrix(m,round(TargetDur(m))) = 1;
    end
end
matrix = conv2(intertapMatrix, convMatrix, 'same');

% Filter
hfilt = ones(window*factor, window)/(factor*window^2);
matrix = imfilter(matrix, hfilt, 'replicate');
    
% Normalize
for m = 1:size(matrix,1)
    matrix(m,:) = matrix(m,:)/sum(matrix(m,:));
end