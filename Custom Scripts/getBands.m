function [feat, featNorm] = getBands(featMat, baseline, starts, ends)
%

recBlock = featMat(baseline, :);
recBlock(recBlock ==1) = NaN;

feat.m(1) = mean(nanmean(recBlock(:, baseline), 1));
feat.std(1) = std(nanmean(recBlock(:, baseline), 1));
for i = 2:length(starts)
    feat.m(i) = mean(nanmean(recBlock(:, starts(i):ends(i)), 1));
    feat.std(i) = std(nanmean(recBlock(:, starts(i):ends(i)), 1));
end

featNorm.m = feat.m./feat.m(1);
featNorm.std = feat.std./feat.m(1);