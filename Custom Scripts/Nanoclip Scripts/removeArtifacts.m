function cleanSnips = removeArtifacts(dirtySnips)
%Stationary, stereotyped artifact removal from stimSnips
cleanSnips = [];

%Cycle through each snip and process
for i = 1:size(dirtySnips,1)
    %Simplify coding
    x = dirtySnips(i,:);
    
    %Find artifact locations via the sharp peaks
    [~, peaksIdx] = findpeaks(x.^2);
    peaksY = x(peaksIdx);
    threshB = peaksY < -5; %Could substitute thresholding on a histogram of amplitudes...
    peaksIdx = peaksIdx(threshB);
    peaksY = peaksY(threshB);
    
%     %Plot to confirm
    figure(999); clf
%     plot(x, 'k'); hold on
%     plot(peaksIdx, peaksY, 'or')
    
    %Estimate artifact spacing
    artSpacing = median(diff(peaksIdx));
    artLength = 0.9 * artSpacing;
    artWing = floor(artLength/2);
    
    %Snip out artifact samples from each instance
    arts = [];
%     figure (9999); clf
    for j = 1:numel(peaksIdx)
        starts = peaksIdx(j) - artWing;
        ends = peaksIdx(j) + artWing;
        
        if starts<1 || ends>numel(x)
            display(['Over ran the file edges for snip ' num2str(j)])
            
        else
            artIdx = starts:ends;
            arts = [arts; x(artIdx)];
%             ts = (1:numel(artIdx))-artWing;
            
%             %Plot to visualize
%             plot(ts, arts(end,:)); hold on
        end
        
    end
    
    %mean artifact waveform
    mArt = mean(arts,1);
%     plot(ts, mArt,'k', 'LineWidth',3)
    
    %Create a filtering/cancelling array
    impArray = zeros(size(x));
    impArray(peaksIdx) = 1;
    filtArray = conv(impArray, -mArt);
    filtArray = filtArray(artWing:end-artWing-1);
    cleanSnips(i,:) = x + filtArray; 
     
%     figure(99999); clf
%     subplot(3,1,1)
%     plot(x)
%     axis tight; ylim([-20, 20])
%     
%     subplot(3,1,2)
%     plot(filtArray)
%     axis tight; ylim([-20, 20])
%     
%     subplot(3,1,3)
%     plot(x+filtArray)
%     axis tight; ylim([-20, 20])
end




