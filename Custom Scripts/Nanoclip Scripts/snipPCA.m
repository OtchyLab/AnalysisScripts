%Code written in response to reviewer request to use more interpretable/repeatable low-D 
%embedding for the fictive singing analysis
%
%Perform pca embedding of spectrogram snips from the fictive singing
%experiment. Script will read in data from a file containing the already
%formatted block of spectrograms, do the pca analysis, plot the
%embedding, and calculate the pairwise distances.
%
%
%Written by TMO 01/27/20

%Clear the workspace
clear all

%Source and destination data
% sLoc = '/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/FigSnips';
sLoc = 'E:\Dropbox\Dropbox'; %for my work desktop

%Load variables from file
load([sLoc, filesep, 'Bird XXX.mat']);
if ~exist('sylType', 'var')
    disp('Uh-oh... something is wrong')
    return
end

%Transform 2-D spectrograms into a 1-D array
specVec = [];
for i = 1:size(sylSpecs,1)
    sp = squeeze(sylSpecs(i, :, :));
    specVec(i,:) = sp(:);
end

%% Perform PCA embedding
rng('default') % for reproducibility

%pca on zscored data
[coeff,score,latent,tsquared,explained,mu] = pca(zscore(specVec));

PC1 = score(:,1);
PC2 = score(:,2);
PC3 = score(:,3);
%% Plot it
figure(45); clf

% Linespec
cols = linspecer(24);
colR = flipud(cols);

un = unique(sylType);
for i = 1:numel(un)
    %Determine which indices are for a given pattern
    idx = find(sylType==un(i));

    %Plot each timeseries in lowD space
%     scatter(Y(idx,1), Y(idx,2), 50, 'MarkerFaceColor', cols(i,:),'MarkerEdgeColor', [0,0,0], 'Marker', 'o'); hold on
    plot(PC1(idx), PC2(idx), 'Marker', 'o', 'Color', cols(i,:), 'LineStyle', 'none', 'MarkerSize', 6); hold on

end
axis square
xlim([-250, 250]); ylim([-200, 200]);
xlabel('PC1'); ylabel('PC2')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [-250,0,250], 'YTick', [-200,0,200])
set(gcf, 'Units', 'Inches', 'Position', [3.5, 3.5, 6, 4.5])


%% Calculate distances in the PC plot above
inMat = []; outMat = [];
for i = 1:numel(un)
    %Determine which indices are for a given pattern
    idx = find(sylType == un(i));
    n_idx = find(sylType ~= un(i));
    
    %In-group distances
    inDist = []; 
    for j = 1:numel(idx)
        U = [PC1(idx(j)), PC2(idx(j)), PC3(idx(j))];
        
        for k = 1:numel(idx)
            if j ~= k
                V = [PC1(idx(k)), PC2(idx(k)), PC3(idx(k))];
                
                D = pdist2(V, U, 'euclidean');
                inDist = [inDist, D];
                
            end
        end
        
    end
    inMat{i} =  inDist;
    
    %Out-group distances
    outDist = [];
    for j = 1:numel(idx)
        U = [PC1(idx(j)), PC2(idx(j)), PC3(idx(j))];
% U = [mean(Y(idx,1)), mean(Y(idx,2))];
        
        for k = 1:numel(n_idx)
            
            V = [PC1(n_idx(k)), PC2(n_idx(k)), PC3(n_idx(k))];
            
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

in = mean(m_in);    %  in = 30.96 (2D) 36.36 (3D)
out = mean(m_out);  % out = 105.72 (2D) 120.87 (3D)




