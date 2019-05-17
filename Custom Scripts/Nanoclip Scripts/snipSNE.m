%Perform t-sne embedding of spectrogram snips from the fictive singing
%experiment. Script will read in data from a file containing the already
%formatted block of spectrograms, do the t-sne analysis, plot the
%embedding, and calculate the pairwise distances.
%
%
%Written by TMO 05/16/19

%Clear the workspace
%clear all

%Source and destination data
% sLoc = '/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/FigSnips';
sLoc = 'E:\Dropbox\Dropbox';

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

%% Perform tSNE embedding
rng('default') % for reproducibility
% options = statset('MaxIter',5000);
[Y, loss] = tsne(specVec, 'Algorithm', 'exact', 'NumPCAComponents', 50, 'Exaggeration', 1, 'Distance', 'seuclidean', 'NumDimensions', 2, 'Perplexity', 35);

%% Plot it
figure(40); clf

% Linespec
cols = linspecer(24);

un = unique(sylType);
for i = 1:numel(un)
    %Determine which indices are for a given pattern
    idx = find(sylType==un(i));

    %Plot each timeseries in lowD space
    scatter(Y(idx,1), Y(idx,2), 50, 'MarkerFaceColor', cols(i,:),'MarkerEdgeColor', [0,0,0], 'Marker', 'o'); hold on

end
axis square
xlim([-50, 50]); ylim([-50, 50]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', [-50,0,50], 'YTick', [-50,0,50])
set(gcf, 'Units', 'Inches', 'Position', [3.5, 3.5, 6, 4.5])


%% Calculate distances in the tSNE plot above
inMat = []; outMat = [];
for i = 1:numel(un)
    %Determine which indices are for a given pattern
    idx = find(sylType == un(i));
    n_idx = find(sylType ~= un(i));
    
    %In-group distances
    inDist = []; 
    for j = 1:numel(idx)
        U = [Y(idx(j),1), Y(idx(j),2)];
        
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

in = mean(m_in);    %  in = 11.6700
out = mean(m_out);  % out = 36.6576




