


% %Inactivation time
% inactStartTime = '11:35';
% sp = regexp(annotName, '_', 'split');
% inactStartDay = sp{2};
% inactDatenum = datenum([inactStartDay '_' inactStartTime], 'yymmdd_HH:MM');
% 
%Indices for selecting pre/post
sylInx = (sylType>=1 & sylType<=100) | sylType==103;
preInx = sylAbsStartTime<inactDatenum;
postInx = sylAbsStartTime>=inactDatenum;

%sub-select the duration and feature arrays
preDur = sylDur(sylInx & preInx);
postDur = sylDur(sylInx & postInx);
preFeat = sylFeat(sylInx & preInx);
postFeat = sylFeat(sylInx & postInx);

%Window Data
winSize = 500;
winStep = 50;

preDurWin = windowTS(preDur, winSize, winStep, 'pad');
preFeatWin = windowTS(preFeat, winSize, winStep, 'pad');

postDurWin = windowTS(postDur, winSize, winStep, 'pad');
postFeatWin = windowTS(postFeat, winSize, winStep, 'pad');

M = []; %Movie file output
m = 1;
preStep = size(preDurWin, 1);
postStep = size(postDurWin, 1);
post1Step = size(post1DurWin, 1);
map = [];
for i = 1:preStep
    figure(20)
    [edgesX, edgesY, N, ~] = ndhist(preFeatWin(i, :), preDurWin(i, :), 'bins', 6, 'filter', 'axis', [-4.5 0 0 450]);
    
    interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
    interpN(isnan(interpN)) = 0;
    norminterpN = interpN./(sum(interpN(:)));
    
    figure(22)
    imagesc(linspace(-4, 0,100), linspace(0, 400, 150), norminterpN, [0, 0.0025]);
    xlabel('Log Entropy')
    ylabel('Syllable Duration (ms)')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14)
    axis xy
    drawnow
    map(m,:,:) = norminterpN;
    
    %axis 'tight'
    if isempty(M)
        M = getframe(gcf);
    else
        M(m) = getframe(gcf);
    end
    m = m + 1;
end


for i = 1:postStep
    figure(20)
    [edgesX, edgesY, N, h] = ndhist(postFeatWin(i, :), postDurWin(i, :), 'bins', 6, 'filter', 'axis', [-4.5 0 0 450]);
    
    interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
    interpN(isnan(interpN)) = 0;
    norminterpN = interpN./(sum(interpN(:)));
    
    figure(22)
    imagesc(linspace(-4, 0,100), linspace(0, 400, 150), norminterpN, [0, 0.0025])
    xlabel('Log Entropy')
    ylabel('Syllable Duration (ms)')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14)
    axis xy
    drawnow
    map(m,:,:) = norminterpN;
    
    %axis 'tight'
    M(m) = getframe(gcf);
    m = m + 1;
end

for i = 1:post1Step
    figure(20)
    [edgesX, edgesY, N, h] = ndhist(post1FeatWin(i, :), post1DurWin(i, :), 'bins', 6, 'filter', 'axis', [-4.5 0 0 450]);
    
    interpN = interp2(edgesX, edgesY, N, linspace(-4, 0,100)', linspace(0, 400, 150));
    interpN(isnan(interpN)) = 0;
    norminterpN = interpN./(sum(interpN(:)));
    
    figure(22)
    imagesc(linspace(-4, 0,100), linspace(0, 400, 150), norminterpN, [0, 0.0025])
    xlabel('Log Entropy')
    ylabel('Syllable Duration (ms)')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 14)
    axis xy
    drawnow
    map(m,:,:) = norminterpN;
    
    %axis 'tight'
    M(m) = getframe(gcf);
    m = m + 1;
end

% movie(M);
movie2avi(M, 'FeatureDriftMoviePre.avi')


for i = 1:size(map,1)
    sq = squeeze(map(i,:,:));
    mapsArray(i,:) = sq(:);
end


euclidDist = squareform(pdist(mapsArray,'euclidean'));
figure
imagesc(euclidDist)
recoveryBlock = euclidDist(1:41, :);
recoveryBlock(recoveryBlock ==0) = NaN;
recovery = nanmean(recoveryBlock,1);
figure
plot(recovery)