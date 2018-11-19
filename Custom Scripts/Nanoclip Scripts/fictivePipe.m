% Script for sequencing the fictive singing analysis
%This operates on the outputs of makeStimSnips.m; principle output (at this
%point) are plots of syllable trajectories through PCA space.

bTSNE = false;

%Load the required snips packages from file

motherLoc = 'C:\Users\Tim\Desktop\Anestetized Singing Bird\Analyzed Data\2018-05-08 LW28\';
toLoad = {'snipsPat11.mat', 'snipsPat12.mat','snipsPat15.mat', 'snipsPat16.mat','snipsPat17.mat', 'snipsPat18.mat','snipsPat19.mat', 'snipsPat20.mat'};
birdMarker = 1;
% motherLoc = 'C:\Users\Tim\Desktop\Anestetized Singing Bird\Analyzed Data\2018-03-30 RY176\';
% toLoad = {'snipsPat01.mat', 'snipsPat03.mat','snipsPat05.mat', 'snipsPat13.mat'};
% birdMarker = 2;

%Define/clear structures
%featCube = []; snipsCube = []; startsCube = []; filesCube = []; patternList = []; snipSet = []; zCube = [];birdID = [];

%Sequentially load files for processing
for i = 1:numel(toLoad)
    %Load from file
    load([motherLoc toLoad{i}]);
    
    numSyls = size(snips,1);
    patNum = toLoad{i}(9:10);
    patName = repmat(str2num(patNum), [numSyls, 1]);
    birdName = repmat(birdMarker, [numSyls, 1]);
    
    %Stack them in a matrix
    snipsCube = [snipsCube; snips];
    startsCube = [startsCube; starts];
    filesCube = [filesCube; files];
    patternList = [patternList; patName];
    birdID = [birdID; birdName];
end

%%
%Calculate feature projections for each snip in the
for i = 1:size(snipsCube,1)
    %extract features
%     S = featureBuilder(snipsCube(i,:));
%     S = S(:,50:250);
    
    %Spec-based prep
    S = specBuilder(snipsCube(i,:));
    padmat = gaussPad(S);
    S = [S;padmat];
    
    %Add to the PCA list
    snipSet = [snipSet, S];
end

%Normalize (z-score) each of the features over all snips and timepoints
snipSet(isnan(snipSet)) = 350;
zSnips = zscore(snipSet');
    
%PCA analysis
[coeff, score, latent, ~, explained] = pca(zSnips);

%Unwrap the scores matrix into syllables again
snipSize = size(score,1)/size(snipsCube,1);
sc = score';
for i = 1:size(snipsCube,1)
    starts = snipSize*(i-1)+1;
    ends = snipSize*(i);
    featCube(i,:,:) = sc(:, starts:ends);
end

%Plotting
t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with
patternTypes = unique(patternList);

%% 
%This is a plot of first 2 or 3 PCs for every syllable token
%Nice plot but a little busy
figure(668); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

    for j = idx'
        %Plot first 3 PCs in 3D space
        PC1 = smooth(squeeze(featCube(j,1,:))', 5, 'moving')';
        PC2 = smooth(squeeze(featCube(j,2,:))', 5, 'moving')';
        PC3 = smooth(squeeze(featCube(j,3,:))', 5, 'moving')';
        plot3(PC1, PC2, PC3, 'Color', colors(i,:), 'Marker', '.', 'LineWidth', 0.25, 'MarkerSize', 6); hold on
%         plot(PC1, PC2, 'Color', colors(i,:)); hold on
    end
end
axis tight
lims = axis; [az, el] = view;
set(gca, 'Box', 'off', 'TickDir', 'out')
%set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3],'ZTick', [-2,0,2])
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
grid
%%

% % %This is a plot of first 3 PCs averages over syllable types
% % %Very clear presentation
% % figure(67); clf
% % for i = 1:numel(patternTypes)
% %     %Determine which indices are for a given pattern
% %     idx = find(patternList==patternTypes(i));
% % 
% %     %Plot first 3 PCs in 3D space
% % %     PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
% % %     PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
% % %     PC3 = smooth(squeeze(mean(featCube(idx,3,:)))', 5, 'moving')';
% %     PC1 = squeeze(mean(featCube(idx,1,:)));
% %     PC2 = squeeze(mean(featCube(idx,2,:)));
% %     PC3 = squeeze(mean(featCube(idx,3,:)));
% %     plot3(PC1, PC2, PC3, 'Color', colors(i,:), 'Marker', '.', 'LineWidth', 1, 'MarkerSize', 10); hold on
% % end
% % % xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
% % axis(lims)
% % set(gca, 'Box', 'off', 'TickDir', 'out')
% % set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3],'ZTick', [-2,0,2])
% % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% 
% 
% % %This is a plot of first 2 PCs averages over syllable types
% % %Very clear presentation
% % figure(670); clf
% % for i = 1:numel(patternTypes)
% %     %Determine which indices are for a given pattern
% %     idx = find(patternList==patternTypes(i));
% % 
% %     %Plot first 3 PCs in 3D space
% %     PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
% %     PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
% %     plot(PC1, PC2, 'Color', colors(i,:), 'Marker', '.', 'LineWidth', 1, 'MarkerSize', 10); hold on
% % end
% % % xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
% % %axis(lims)
% % set(gca, 'Box', 'off', 'TickDir', 'out')
% % set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3],'ZTick', [-2,0,2])
% % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% 
% % % %This is a plot of first 2 PCs averages over syllable types
% % % %Not correct and very misleading
% % % figure(672); clf
% % % for i = 1:numel(patternTypes)
% % %     %Determine which indices are for a given pattern
% % %     idx = find(patternList==patternTypes(i));
% % % 
% % %     %Plot first 3 PCs in 3D space
% % %     PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
% % %     PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
% % %     PC1s = smooth(squeeze(std(featCube(idx,1,:))./sqrt(numel(idx)))', 5, 'moving')';
% % %     PC2s = smooth(squeeze(std(featCube(idx,2,:))./sqrt(numel(idx)))', 5, 'moving')';
% % %     
% % %     plot(PC1, PC2, 'Color', colors(i,:), 'Marker', '.', 'LineWidth', 1, 'MarkerSize', 10); hold on
% % %     %plot(PC1+PC1s, PC2+PC2s, 'Color', colors(i,:), 'LineWidth', 1);
% % %     %plot(PC1-PC1s, PC2-PC2s, 'Color', colors(i,:), 'LineWidth', 1);
% % %     patch([PC1+PC1s, fliplr(PC1-PC1s)],[PC2+PC2s, fliplr(PC2-PC2s)], colors(i,:), 'FaceAlpha', .3, 'EdgeColor','none')
% % % end
% % % % xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
% % % axis(lims)
% % % set(gca, 'Box', 'off', 'TickDir', 'out')
% % % set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3])
% % % xlabel('PC1'); ylabel('PC2');
% % 
% % 
% % % %This is a plot of first 3 PCs averages over syllable types, with plane for
% % % %sem -- not a great image.
% % % figure(68); clf
% % % for i = 1:numel(patternTypes)
% % %     %Determine which indices are for a given pattern
% % %     idx = find(patternList==patternTypes(i));
% % % 
% % %     %Plot first 3 PCs in 3D space
% % %     PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
% % %     PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
% % %     PC3 = smooth(squeeze(mean(featCube(idx,3,:)))', 5, 'moving')';
% % %     PC1s = smooth(squeeze(std(featCube(idx,1,:))./numel(idx))', 5, 'moving')';
% % %     PC2s = smooth(squeeze(std(featCube(idx,2,:))./numel(idx))', 5, 'moving')';
% % %     PC3s = smooth(squeeze(std(featCube(idx,3,:))./numel(idx))', 5, 'moving')';
% % % %     PC1 = squeeze(mean(featCube(idx,1,:)));
% % % %     PC2 = squeeze(mean(featCube(idx,2,:)));
% % % %     PC3 = squeeze(mean(featCube(idx,3,:)));
% % %     plot3(PC1, PC2, PC3, 'Color', colors(i,:), 'Marker', '.', 'LineWidth', 1, 'MarkerSize', 10); hold on
% % %     fill3([PC1+PC1s, fliplr(PC1-PC1s)],[PC2+PC2s, fliplr(PC2-PC2s)],[PC3+PC3s, fliplr(PC3-PC3s)], colors(i,:), 'FaceAlpha', .3, 'EdgeColor','none')
% % % end
% % % % xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
% % % axis(lims)
% % % set(gca, 'Box', 'off', 'TickDir', 'out')
% % % set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3],'ZTick', [-2,0,2])
% % % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% % % %axis(l); view([az,el])
% % 
% % 
% % %This is a plot of first 3 PCs averages over syllable types
% % %with circles indicating the mean euclidean distance at each 1ms time point
% % figure(90); clf
% % PCDist = [];
% % for i = 1:numel(patternTypes)
% %     %Determine which indices are for a given pattern
% %     idx = find(patternList==patternTypes(i));
% % 
% %     %Plot first 3 PCs in 3D space
% %     PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
% %     PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
% %     PC3 = smooth(squeeze(mean(featCube(idx,3,:)))', 5, 'moving')';
% %     V = [PC1', PC2', PC3'];
% %     Ds = [];
% %     
% %     %Calculate the point-wise Euclidean distance for each syllable
% %     for j = idx'
% %         PC1i = smooth(squeeze(featCube(j,1,:))', 5, 'moving')';
% %         PC2i = smooth(squeeze(featCube(j,2,:))', 5, 'moving')';
% %         PC3i = smooth(squeeze(featCube(j,3,:))', 5, 'moving')';
% % 
% %         %Pointwise EU Dist
% %         D = pdist2(V, [PC1i', PC2i', PC3i'], 'euclidean');
% %         Ds(j,:) = diag(D)';
% %         
% %     end
% %     
% %     %Calc mean pointwise distance
% %     PCDist(i,:) = mean(Ds,1);
% %     PCDistRad = (PCDist(i,:).^2).*pi; %Radius to area conversion
% %     
% %     %Plot
% %     plot3(PC1, PC2, PC3, 'Color', colors(i,:), 'LineWidth', 1); hold on
% %     scatter3(PC1, PC2, PC3, 100*PCDistRad, colors(i,:))
% %     
% % end
% % % xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
% % axis(lims)
% % set(gca, 'Box', 'off', 'TickDir', 'out')
% % set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3],'ZTick', [-2,0,2])
% % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% %

%%
%This is a plot of first 2 PCs averages over syllable types
%Shaded band indicating mean euclidean distance
figure(71); clf
PCDist = [];
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

    %Plot first 3 PCs in 3D space
    PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
    PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';

    V = [PC1', PC2'];
    Ds = [];
    
    %Calculate the point-wise Euclidean distance for each syllable
    for j = idx'
%         PC1i = smooth(squeeze(featCube(j,1,:))', 5, 'moving');
%         PC2i = smooth(squeeze(featCube(j,2,:))', 5, 'moving');
        PC1i = squeeze(featCube(j,1,:));
        PC2i = squeeze(featCube(j,2,:));
        
        %Pointwise EU Dist
        D = pdist2(V, [PC1i, PC2i], 'euclidean');
        Ds(j,:) = diag(D)';
    end
    
    %Calc mean pointwise distance
    PCDist(i,:) = mean(Ds,1)./5;
    [innerx, innery, outerx, outery] = xyshadedbar(PC1, PC2, PCDist(i,:));
    
    %Plot
    plot(PC1, PC2, 'Color', colors(i,:), 'LineWidth', 2); hold on
    
    %Lines
    plot(innerx, innery, 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 0.5);
    plot(outerx, outery, 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 0.5);
    
    %Shaded
%     patch([innerx, fliplr(outerx)], [innery, fliplr(outery)], colors(i,:), 'FaceAlpha', .3, 'EdgeColor','none')

end
% xlim([-4, 4]); ylim([-3, 3]); zlim([-2.5, 2.5]);
set(gca, 'Box', 'off', 'TickDir', 'out')
%set(gca, 'XTick', [-3,0,3],'YTick', [-3,0,3])
xlabel('PC1'); ylabel('PC2');


%% 
%This is a plot of first 3 PCs averages vs time  over syllable types
%Shaded band indicating mean euclidean distance
figure(78); clf
PCDist = [];
% timepnts = -100:392;
timepnts = (1:size(featCube,3))-40;
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

    %Plot first 3 PCs in 3D space
    PC1 = smooth(squeeze(mean(featCube(idx,1,:)))', 5, 'moving')';
    PC2 = smooth(squeeze(mean(featCube(idx,2,:)))', 5, 'moving')';
    PC3 = smooth(squeeze(mean(featCube(idx,3,:)))', 5, 'moving')';
    PC1s = smooth(squeeze(std(featCube(idx,1,:), 1, 1))./sqrt(numel(idx))', 5, 'moving')';
    PC2s = smooth(squeeze(std(featCube(idx,2,:), 1, 1))./sqrt(numel(idx))', 5, 'moving')';
    PC3s = smooth(squeeze(std(featCube(idx,3,:), 1, 1))./sqrt(numel(idx))', 5, 'moving')';
    
    V = [PC1, PC2, PC3];
    Ds = [];
    
%     %Calculate the point-wise Euclidean distance for each syllable
%     for j = idx'
%         %         PC1i = smooth(squeeze(featCube(j,1,:))', 5, 'moving');
%         %         PC2i = smooth(squeeze(featCube(j,2,:))', 5, 'moving');
%         PC1i = squeeze(featCube(j,1,:));
%         PC2i = squeeze(featCube(j,2,:));
%         PC3i = squeeze(featCube(j,3,:));
%         
%         %Pointwise EU Dist
%         D = pdist2(V, [PC1i, PC2i, PC3i], 'euclidean');
%         Ds(j,:) = diag(D)';
%     end
%     
%     %Calc mean pointwise distance
%     PCDist(i,:) = mean(Ds,1);
    
    %Plot
    subplot(3,1,1)
            s = shadedErrorBar(timepnts,PC1, PC1s, 'k', 1); hold on
        s.mainLine.Color = colors(i,:);
        s.edge(1).Color = colors(i,:);
        s.edge(2).Color = colors(i,:);
        s.patch.FaceColor = colors(i,:);
    xlim([-5, 110]);
    set(gca, 'Box', 'off', 'TickDir', 'out')
    ylabel('PC1');
    
    subplot(3,1,2)
                s = shadedErrorBar(timepnts,PC2, PC2s, 'k', 1); hold on
        s.mainLine.Color = colors(i,:);
        s.edge(1).Color = colors(i,:);
        s.edge(2).Color = colors(i,:);
        s.patch.FaceColor = colors(i,:);
    xlim([-5, 110]);
    set(gca, 'Box', 'off', 'TickDir', 'out')
    ylabel('PC2');
    
    subplot(3,1,3)
                s = shadedErrorBar(timepnts,PC3, PC3s, 'k', 1); hold on
        s.mainLine.Color = colors(i,:);
        s.edge(1).Color = colors(i,:);
        s.edge(2).Color = colors(i,:);
        s.patch.FaceColor = colors(i,:);
    xlim([-5, 110]);
    set(gca, 'Box', 'off', 'TickDir', 'out')
    xlabel('Time (ms)'); ylabel('PC3');
end
subplot(3,1,1)
ys = ylim;
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
xlim([-5, 110]);
ylabel('PC1')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [0, 50, 100], 'XTickLabel', [])

subplot(3,1,2)
ys = ylim;
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
xlim([-5, 110]);
ylabel('PC2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [0, 50, 100], 'XTickLabel', [])

subplot(3,1,3)
ys = ylim;
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
xlim([-5, 110]);
ylabel('PC3'); xlabel('Time (ms)')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [0, 50, 100])



%%
%tSNE embedding
if ~bTSNE
    return %Abort if you don't want to run this section
end

%%
snipSet = [];
%Calculate feature projections for each snip in the 
for i = 1:size(snipsCube,1)
    %extract features
    S = specBuilder(snipsCube(i,:));
    %S = featureBuilder(snipsCube(i,:));
%     S = S(:,50:250);
     %S = S(:,100:200);
    padmat = gaussPad(S);
    S = [S;padmat];

    %Add to the PCA list
    snipSet = [snipSet, S];
end

%Perform tSNE embedding
rng('default') % for reproducibility
Y = tsne(snipSet','Algorithm','barneshut', 'Distance', 'euclidean', 'NumDimensions', 2, 'Perplexity', 50);
% Y = tsne(snipSet','Algorithm','barneshut', 'Distance', 'euclidean', 'NumDimensions', 2, 'Perplexity', 100);
% figure; scatter(Y(:,1),Y(:,2))

%Unwrap the embedding matrix into syllable snips again
snipSize = size(Y,1)/size(snipsCube,1);
sc = Y';
featCube = [];
for i = 1:size(snipsCube,1)
    starts = snipSize*(i-1)+1;
    ends = snipSize*(i);
    featCube(i,:,:) = sc(:, starts:ends);
end

%Plotting
t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with
patternTypes = unique(patternList);


%%
%This is a plot
figure(39); clf
scatter(Y(:,1), Y(:,2), '.k')
%xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])

% y = featCube(idx,1,:); y = y(:);
% z = featCube(idx,2,:); z = z(:);
% scatter(y,z, '.b')
%%
%This is a plot
figure(40); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

     %Plot each timeseries in lowD space
    for j = idx'
        plot(squeeze(featCube(j,1,:)), squeeze(featCube(j,2,:)), 'LineStyle', 'none', 'Marker', '.', 'Color', colors(i,:)); hold on
    end
end
%xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])

%% This is a plot
figure(41); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

     %Plot each timeseries in lowD space
    for j = idx'
        d1 = smooth(squeeze(featCube(j,1,:))', 9, 'moving')';
        d2 = smooth(squeeze(featCube(j,2,:))', 9, 'moving')';
        
        plot(d1, d2, 'Color', colors(i,:)); hold on
    end
end
xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])

%% This is a plot
figure(42); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));
    
        d1 = smooth(mean(squeeze(featCube(idx,1,:)),1)', 5, 'moving')';
        d2 = smooth(mean(squeeze(featCube(idx,2,:)),1)', 5, 'moving')';
%     d1 = mean(squeeze(featCube(idx,1,:)),1);
%     d2 = mean(squeeze(featCube(idx,2,:)),1);
    plot(d1, d2, 'Marker', '.', 'Color', colors(i,:), 'LineWidth', 1.5); hold on
end
xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])
%% This is a plot
figure(43); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

     %Plot each timeseries in lowD space
    for j = idx'
        d1 = smooth(squeeze(featCube(j,1,:))', 9, 'moving')';
        d2 = smooth(squeeze(featCube(j,2,:))', 9, 'moving')';
        
        subplot(2,1,1)
        plot(d1,'Color', colors(i,:)); hold on
        
        subplot(2,1,2)
        plot(d2,'Color', colors(i,:)); hold on
    end
end

%%
%This is a plot
figure(44); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));
    
        d1 = smooth(mean(squeeze(featCube(idx,1,:)),1)', 5, 'moving')';
        d2 = smooth(mean(squeeze(featCube(idx,2,:)),1)', 5, 'moving')';
        s1 = smooth(std(squeeze(featCube(idx,1,:)),1,1)'./sqrt(numel(idx)), 5, 'moving')';
        s2 = smooth(std(squeeze(featCube(idx,2,:)),1,1)'./sqrt(numel(idx)), 5, 'moving')';
        
        subplot(2,1,1)
        s = shadedErrorBar(-5:(numel(d1)-6),d1, s1, 'k', 1); hold on
        s.mainLine.Color = colors(i,:);
        s.edge(1).Color = colors(i,:);
        s.edge(2).Color = colors(i,:);
        s.patch.FaceColor = colors(i,:);
        
        subplot(2,1,2)
        s = shadedErrorBar(-5:(numel(d1)-6),d2, s2, 'k', 1); hold on
        s.mainLine.Color = colors(i,:);
        s.edge(1).Color = colors(i,:);
        s.edge(2).Color = colors(i,:);
        s.patch.FaceColor = colors(i,:);
end
subplot(2,1,1)
ys = ylim;
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
xlim([-5, 110]);
ylabel('t-SNE D1')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [0, 50, 100], 'XTickLabel', [])

subplot(2,1,2)
ys = ylim;
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
xlim([-5, 110]);
ylabel('t-SNE D2'); xlabel('time (ms)')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [0, 50, 100])


%% Calculate distances in the tSNE plot above
inMat = []; outMat = [];
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));
    n_idx = find(patternList~=patternTypes(i));
    
    %In-group distances
    inDist = []; 
    for j = 1:numel(idx)
        d1u = smooth(squeeze(featCube(idx(j),1,:))', 5, 'moving');
        d2u = smooth(squeeze(featCube(idx(j),2,:))', 5, 'moving');
        U = [d1u, d2u];
        
        for k = 1:numel(idx)
            if j ~= k
                d1v = smooth(squeeze(featCube(idx(k),1,:))', 5, 'moving');
                d2v = smooth(squeeze(featCube(idx(k),2,:))', 5, 'moving');
                V = [d1v, d2v];
                
                D = pdist2(V, U, 'euclidean');
                inDist = [inDist, mean(diag(D))];
                
            end
        end
        
    end
    inMat{i} =  inDist;
    
    %Out-group distances
    outDist = [];
    for j = 1:numel(idx)
        d1u = smooth(squeeze(featCube(idx(j),1,:))', 5, 'moving');
        d2u = smooth(squeeze(featCube(idx(j),2,:))', 5, 'moving');
        U = [d1u, d2u];
        
        for k = 1:numel(n_idx)
            
            d1v = smooth(squeeze(featCube(n_idx(k),1,:))', 5, 'moving');
            d2v = smooth(squeeze(featCube(n_idx(k),2,:))', 5, 'moving');
            V = [d1v, d2v];
            
            D = pdist2(V, U, 'euclidean');
            outDist = [outDist, mean(diag(D))];
            
        end
        
    end
    outMat{i} =  outDist;
end

%Do stats on the distances
m_in = cellfun(@mean, inMat);
s_in = cellfun(@std, inMat);
m_out = cellfun(@mean, outMat);
s_out = cellfun(@std, outMat);

figure(6); clf
xs = 1:2;
for i = 1:numel(m_in)
    %Generate a lkittle noise to jigger the points
   jig = randn(1,2)./25;
   
   %Plot mean and std for each stimulation pattern
   errorbar(xs+jig, [m_in(i), m_out(i)], [s_in(i), s_out(i)], 'Color', colors(i,:), 'Marker', 'o', 'LineStyle', '-'); hold on
    
    
end

%Plot summary stats
bar(xs, [mean(m_in), mean(m_out)], 'k', 'FaceAlpha', 0.3)
errorbar(xs,[mean(m_in), mean(m_out)],[std(m_in), std(m_out)], '.k', 'LineWidth', 2.5)
xlim([0, 3])
ylabel('Mean Pairwise Distance')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Same'; 'Different'})
set(gcf, 'Units', 'Inches', 'Position', [2.7500, 7, 3.5, 5.5])


%% tSNE with the whole spectrogram timeseries condensed to a single point
snipSet = [];
%Calculate feature projections for each snip in the 
for i = 1:size(snipsCube,1)
    %extract features
    S = specBuilder(snipsCube(i,:));
    padmat = gaussPad(S);
    S = [S;padmat];

    %Add to the PCA list
    snipSet = [snipSet, S(:)];
end

%Perform tSNE embedding
rng('default') % for reproducibility
% options = statset('MaxIter',5000);
[Y, loss] = tsne(snipSet','Algorithm','exact', 'NumPCAComponents', 50, 'Distance', 'euclidean', 'NumDimensions', 2, 'Perplexity', 15);
% Y = tsne(snipSet','Algorithm','barneshut', 'Distance', 'euclidean', 'NumDimensions', 2, 'Perplexity', 100);
% figure; scatter(Y(:,1),Y(:,2))

% %Unwrap the embedding matrix into syllable snips again
% snipSize = size(Y,1)/size(snipsCube,1);
% sc = Y';
% featCube = [];
% for i = 1:size(snipsCube,1)
%     starts = snipSize*(i-1)+1;
%     ends = snipSize*(i);
%     featCube(i,:,:) = sc(:, starts:ends);
% end

%Plotting
t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with
patternTypes = unique(patternList);

%%
%This is a plot
figure(39); clf
scatter(Y(:,1), Y(:,2), '.k')
%xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
%set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])

% y = featCube(idx,1,:); y = y(:);
% z = featCube(idx,2,:); z = z(:);
% scatter(y,z, '.b')
%%
%This is a plot
figure(40); clf
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));

     %Plot each timeseries in lowD space
    scatter(Y(idx,1), Y(idx,2), 'MarkerFaceColor', colors(i,:),'MarkerEdgeColor', colors(i,:), 'Marker', 'o'); hold on

end
xlim([-50, 50]); ylim([-50, 50]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-50, 0, 50], 'YTick', [-50, 0, 50])

%% 
figure(939); clf
gscatter(Y(:,1), Y(:,2), patternList)
%xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')
%set(gca, 'XTick', [-60, 0, 60], 'YTick', [-80, 0, 60])

% y = featCube(idx,1,:); y = y(:);
% z = featCube(idx,2,:); z = z(:);
% scatter(y,z, '.b')

%% Calculate distances in the tSNE plot above
inMat = []; outMat = [];
for i = 1:numel(patternTypes)
    %Determine which indices are for a given pattern
    idx = find(patternList==patternTypes(i));
    n_idx = find(patternList~=patternTypes(i));
    
    %In-group distances
    inDist = []; 
    for j = 1:numel(idx)

        U = [Y(idx(j),1), Y(idx(j),2)];
%         U = [mean(Y(idx,1)), mean(Y(idx,2))];
        
        for k = 1:numel(idx)
            if j ~= k
                V = [Y(idx(k),1), Y(idx(k),2)];
                
                D = pdist2(V, U, 'euclidean');
                inDist = [inDist, D];
                
            end
        end
        
    end
    inMat{i} =  inDist;
    
    %Out-group distances
    outDist = [];
    for j = 1:numel(idx)
        U = [Y(idx(j),1), Y(idx(j),2)];
% U = [mean(Y(idx,1)), mean(Y(idx,2))];
        
        for k = 1:numel(n_idx)
            
            V = [Y(n_idx(k),1), Y(n_idx(k),2)];
            
            D = pdist2(V, U, 'euclidean');
            outDist = [outDist, D];
            
        end
        
    end
    outMat{i} =  outDist;
end

%Do stats on the distances
m_in = cellfun(@mean, inMat);
s_in = cellfun(@std, inMat);
m_out = cellfun(@mean, outMat);
s_out = cellfun(@std, outMat);

figure(6); clf
xs = 1:2;
for i = 1:numel(m_in)
    %Generate a lkittle noise to jigger the points
   jig = randn(1,2)./25;
   
   %Plot mean and std for each stimulation pattern
   errorbar(xs+jig, [m_in(i), m_out(i)], [s_in(i), s_out(i)], 'Color', colors(i,:), 'Marker', 'o', 'LineStyle', '-'); hold on
    
    
end

%Plot summary stats
bar(xs, [mean(m_in), mean(m_out)], 'k', 'FaceAlpha', 0.3)
errorbar(xs,[mean(m_in), mean(m_out)],[std(m_in), std(m_out)], '.k', 'LineWidth', 2.5)
xlim([0, 3])
ylabel('Mean Pairwise Distance')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Same'; 'Different'})
set(gcf, 'Units', 'Inches', 'Position', [2.7500, 7, 3.5, 5.5])
















