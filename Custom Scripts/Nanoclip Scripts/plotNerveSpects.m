%Script plots examples of the pre- and post-implant/manipulation song
%spectrograms for the the revised nanoclip paper. As of the writing of this
%description, it will generate the spectrograms that appear in Sup Fig S5.
%
% The script reads in the rawdata from a previously aligned Talon datasets containing 
% songs from the pre and manip + 1d .mat files saved in the spot workspace
% share.
%
% Written by TMO 02-20-20

%Clear it all
clear all

%File locations and names
mother = 'V:\SongbirdData\Analysis';
folder = 'llb59rblk132'; manip = 'clip';
% folder = 'lr11ry176'; manip = 'clip';
% folder = 'lr14rblk17'; manip = 'clip';
% folder = 'LLY12'; manip = 'sham';
% folder = 'LLY13'; manip = 'sham';
% folder = 'LLY14'; manip = 'crush';
% folder = 'LLY16'; manip = 'sham';
% folder = 'LLR04'; manip = 'crush';
% folder = 'LLY72'; manip = 'crush';
% folder = 'LW58'; manip = 'intact';
% folder = 'LW60'; manip = 'intact';
%  folder = 'LY80'; manip = 'intact';
sourceDir = [mother, filesep, folder];

%Assemble file paths
datafiles = dir([sourceDir, filesep, '*dataset.mat']);
files = getStructField(datafiles, 'name');
files = sort(files);

if numel(files) > 3
    display('Hmmm, unexpected number of dataset files here. Double check that this is all ok.')
end

%Set up labels and pointers
pntr = 21;
labels = {'Pre', 'Post'};
t = ['Spectrograms for ' folder ' and pntr=' num2str(pntr)];

%Generate the figures
h = figure(100); clf
set(gcf, 'Units', 'inches', 'Position', [25, 7.5, 5, 4.5])

% Pre spec
fullPath = [sourceDir, filesep, files{1}];
S = load(fullPath, 'audio');

subplot(2,1,1)
imagesc(-squeeze(S.audio.aligned_audioCube(pntr,:,:)), [-6, -1]); colormap(jet)
axis xy; axis tight

set(gca, 'Box', 'off', 'TickDir', 'out', 'YTickLabels', '', 'LineWidth', 1)
xlabel('Time (ms)'); ylabel(labels{1}); title(t)

clear('S');

% Post spec
fullPath = [sourceDir, filesep, files{1}];
S = load(fullPath, 'audio');

subplot(2,1,2)
imagesc(-squeeze(S.audio.aligned_audioCube(pntr,:,:)), [-6, -1]); colormap(jet)
axis xy; axis tight

set(gca, 'Box', 'off', 'TickDir', 'out', 'YTickLabels', '', 'LineWidth', 1)
xlabel('Time (ms)'); ylabel(labels{2})

clear('S');

%Save the plot to file
savefig(h, [sourceDir, filesep, 'prepostSpecs pntr=' num2str(pntr) '.fig'])
