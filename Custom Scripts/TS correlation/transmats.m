function out = transmats(audio,neuro,template)

out.NSmAP = mean(audio);
out.NSsAP = std(audio);

out.NSmNP = mean(neuro);
out.NSsNP = std(neuro);

[Dval out.p] = DTWsubWeightedNormBand(template,mean(audio));
out.NSmAP = alignSeries(mean(audio),p);
out.NSsAP = alignSeries(std(audio),p);

out.mNP = alignRawTS(out.NSmNP,length(mean(audio)),p);
out.sNP = alignRawTS(out.NSsNP,length(mean(audio)),p);

% out.NSmAP = mean(audio);
% out.NSsAP = std(audio);
% 
% out.mNP = mean(neuro);
% out.sNP = std(neuro);

target = length(out.mNP);
out.mAP = interp1(1:target/length(audio):target,mean(audio),1:target,'linear');
out.sAP = interp1(1:target/length(audio):target,std(audio),1:target,'linear');

