%Perform t-sne embedding of spectrogram snips from the fictive singing
%experiment. Script will read in data from a file containing the already
%formatted block of spectrograms, do the t-sne analysis, plot the
%embedding, and calculate the pairwise distances.
%
%
%Written by TMO 05/16/19

%Clear the workspace
clear all

%Source and destination data
sLoc = '/Users/Tim/Desktop/Nanoclip Paper Figure Stuff/FigSnips';

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
[Y, loss] = tsne(specVec','Algorithm','exact', 'NumPCAComponents', 50, 'Distance', 'euclidean', 'NumDimensions', 2, 'Perplexity', 15);


%Plot it?
figure(939); clf
gscatter(Y(:,1), Y(:,2), patternList)
%xlim([-60, 60]); ylim([-80, 60]);
xlabel('t-SNE D1'); ylabel('t-SNE D2')
set(gca, 'Box', 'off', 'TickDir', 'out')








