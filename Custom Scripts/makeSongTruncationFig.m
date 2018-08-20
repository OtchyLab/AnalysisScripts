%Test script to create the song truncation figures for the Inactivations paper

%clear all

%Define the source material
annotName = 'Grn046_140518_annotation.mat';           %Annotation to process
critTime = datenum(2014, 5, 15, 11, 32, 00);                %Lesion time and date
motif = [1,2,3];                                                           %sequence defined in the annotation

%Clear/declare variables
fs = 44150;                 %in Hz
buffer = 0.010;             %in seconds
maxSnip = -inf;          %initialize to -inf
outputSnips = [];
outputSylNums = [];
outputPos = [];
outputFile = [];
outputLength = [];
outputTime = [];

%Generate all other file locations/names based on annotation name
mother = 'V:\';
sp = regexp(annotName, '_', 'split');
annotLoc = [mother, char(sp{1}) filesep annotName];
fileSource = [mother, char(sp{1}) filesep '20' char(sp{2}(1:2)) '-' char(sp{2}(3:4)) '-' char(sp{2}(5:6)) filesep];

%Load the annotation from file
load(annotLoc)

%Step through each recording file and analyze independently
% for i = 1:length(elements)
for i = 1:length(elements)
    %Extract from the annotation for clean coding
    starts = elements{i}.segFileStartTimes;
    ends = elements{i}.segFileEndTimes;
    types = elements{i}.segType';
    recDateNum = getFileTime(keys{i});
    
    %Find sequences in the current file that correspond to the motif or subsets thereof
    motifRaw = [];
    for j = 1:length(motif)
        motifRaw{j} = strfind(types,motif(1:j));
    end
    
    %Parse these sequences for uniqueness -- longer motif gets to dominate
    motifTrim{j} = motifRaw{j};
    for j = length(motif):-1:2
        %Locate and remove common starts for shorter sequences
        overlap = ~ismember(motifRaw{j-1}, motifRaw{j});
        motifTrim{j-1} = motifRaw{j-1}(overlap);
    end

    %Assembe the song snip times
    snipTimes = [];
    snipPos = [];
    snipLength = [];
    for j = 1:length(motif)
        if ~isempty(motifTrim{j})
            for k = 1:length(motifTrim{j})
                idx = motifTrim{j}(k):(motifTrim{j}(k) + j-1);
                snipTimes = [snipTimes; starts(idx(1))-buffer, ends(idx(end))+buffer];
                snipPos = [snipPos; motifTrim{j}(k)];
                snipLength = [snipLength; j];
            end
        end
    end
    
    %Sort the snip times in the order they appear in the recording
    [snipPos, I] = sort(snipPos);
    snipTimes = snipTimes(I,:);
    snipLength = snipLength(I);
    
    %If there are snips to make, load the file into memory and extract the snip
    if ~isempty(snipTimes)
        %Load audio
        if strcmp(keys{i}(end), 'v')
            %It's a WAV
            audio = audioread([fileSource keys{i}]);
        else
            %It's a DAT
            [chans, ~] = getChannels([fileSource keys{i}]);
            audio = chans(1,:);
            clear('chans');
        end
        
        %Snip extraction
        for j = 1:length(snipPos)
            snipIdx = floor(snipTimes(j,1)*fs):ceil(snipTimes(j,2)*fs);
            snip = audio(snipIdx);
            maxSnip = max([maxSnip,length(snip)]);
            
            %Shuffle the snip and associated data off to output arrays
            outputSnips{end+1} = snip;
            outputLength(end+1) = snipTimes(j,2) - snipTimes(j,1);
            outputSylNums(end+1) = snipLength(j);
            outputPos(end+1) = snipPos(j);
            outputFile{end+1} = keys{i};
            outputTime (end+1) = recDateNum;
        end
    end
    
end

%Now that the snips are assembled, they need to be displayed in some way

%Use a flattened spectrogram (5ms running window) to smooth and flatten audio
[~, ~, ~, p] = cellfun(@(x) spectrogram((x/(sqrt(mean(x.^2)))),220,220-44,512,44150), outputSnips, 'UniformOutput', 0);
q = cellfun(@(x) mat2gray(sum(log10(abs(x)),1)), p, 'UniformOutput', 0);


%Calculate threshold for alignment
% maxSpec = max(cell2mat(cellfun(@(x) length(x), q, 'UniformOutput', 0)));
maxSpec = 550;
s = cellfun(@(x) [x, zeros(1,maxSpec-length(x))], q, 'UniformOutput', 0)';
thresh = 0.35*max(mean(cell2mat(s),1));

%use threshold-crossing to align to the onsetr of the first syllable and make each the same length by padding with zeros
alignedSpecs = cellfun(@(x) [zeros(1,100-find(x>thresh,1,'first')), x, zeros(1,maxSpec+find(x>thresh,1,'first')-length(x))], q, 'UniformOutput', 0);
finalLength = cell2mat(cellfun(@(x) length(x)-find(x>thresh,1,'first')-(length(x)-find(x>thresh,1,'last')), q, 'UniformOutput', 0));

%And convert to a uniform matrix
SpecsMat = cell2mat(alignedSpecs');



index = sort(randi([400,1200],300,1)); %For the first day
% index = sort(randi(length(finalLength),300,1));
% index = 50:351; index(index==86) = [];

%Plot the aligned amplitude of the motifs
% ticks = (1:100:size(SpecsMat,2)) + 100;
% tickLab = 0:100:(size(SpecsMat,2));
% h(1) = figure(100); clf
% imagesc(SpecsMat); colormap(jet)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ticks, 'XTickLabel', tickLab)
% xlabel('Time (ms)')
% ylabel('Renditions')
% title('Motif Amplitudes')

%Bar plot of segment lengths
h(2) = figure(202); clf
for k = 1:length(index)%length(finalLength)
    b(k) = barh(k, finalLength(index(k)), 1, 'EdgeColor', 'none'); hold on
    if outputSylNums(index(k)) == 3
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,400])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:100:600);
xlabel('Song Length (ms)')
ylabel('Renditions')


%Select a subset for the figure
% subSpec.D2 = SpecsMat(index,:);
% subLength.D2 = finalLength(index);

%%
%For plotting the summary figure
figure(200); clf
subplot(5,1,1)
imagesc(subSpec.pre); cs = caxis;
title('Motif Amplitudes')
ylabel('Pre Lesion')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(5,1,2)
imagesc(subSpec.post, cs);
ylabel('Post Lesion')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(5,1,3)
imagesc(subSpec.PD1, cs);
ylabel('Post D1')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(5,1,4)
imagesc(subSpec.PD2, cs);
ylabel('Post D2')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(5,1,5)
imagesc(subSpec.PD3, cs); colormap(jet)
ticks = (1:100:size(SpecsMat,2)) + 100;
tickLab = 0:100:(size(SpecsMat,2));
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', ticks, 'XTickLabel', tickLab)
xlabel('Time (ms)')
ylabel('Post D3')


figure(201); clf
subplot(4,1,1)
for k = 1:length(subLength.pre)
    b(k) = barh(k, subLength.pre(k), 1, 'EdgeColor', 'none'); hold on
    if subLength.pre(k) > 275
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,400])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:100:600);

subplot(4,1,2)
for k = 1:length(subLength.post)
    b(k) = barh(k, subLength.post(k), 1, 'EdgeColor', 'none'); hold on
    if subLength.post(k) > 275
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,400])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:100:600);

subplot(4,1,3)
for k = 1:length(subLength.D1)
    b(k) = barh(k, subLength.D1(k), 1, 'EdgeColor', 'none'); hold on
    if subLength.D1(k) > 275
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,400])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:100:600);

subplot(4,1,4)
for k = 1:length(subLength.D2)
    b(k) = barh(k, subLength.D2(k), 1, 'EdgeColor', 'none'); hold on
    if subLength.D2(k) > 275
        b(k).FaceColor = [0.5, 0.5, 0.5];
    else
        b(k).FaceColor = 'r';
    end
end
axis tight; axis ij
xlim([0,400])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:100:600);
