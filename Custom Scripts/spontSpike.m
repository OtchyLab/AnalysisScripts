function presentCell = spontSpike

%%%%%%%%%%%%%%%%%%%%%%%
% File setup
%%%%%%%%%%%%%%%%%%%%%%%

%Set DAT source directory and index DAT files
sourceDir = 'V:\Grn141\Spon\2014-12-10';
filelist = dir([sourceDir filesep '*.dat']);
for i = 1:length(filelist)
    %Extract out the filename
    fileNames(i,:) = filelist(i).name;
end

%Load (specific) annotation record
annotName = 'V:\Grn141\Grn141_141210_annotationSPON.mat';
load(annotName,'elements', 'keys')
elements = elements;
keys = keys;

%Set some constants
fs = 44150;
chan = 4;
n_std = 4; %number of standard deviations for detection threshold
direction = 'pos'; %apply the spike threshold to the positive/negative side

%Set filtering params
HP_fNorm = 1000/(fs/2); % 1kHz
LP_fNorm = 15000/(fs/2); % 15kHz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]); %Butterworth, 2-p

%Initialize arrays
numRec = length(elements);
presentCell = [];

%Cycle through each record in the annotation and process the corresponding DAT-file.
for i = 1:numRec
    if ~ismember(keys{i},fileNames)
        %If it's not in the source dir, throw and error message
        display([keys{i} ' not found in the source folder'])
    else
        %If it is, load the file into the workspace, parse for the selected channel, and filter as specified
        fullName = [sourceDir filesep keys{i}];
        [rawdata, ~] = getChannels(fullName);
        data = filtfilt(BP_b,BP_a, rawdata(chan+1, :));
        
        %And get the file timestamp while you're here
        fileTimes(i) = getFileTime(keys{i});
    
        %Calculate threshold for detection
        d_m = mean(data);
        d_std = std(data);
        if strcmp(direction, 'pos')
            thresh(i) = d_m + (n_std*d_std);
        elseif strcmp(direction, 'neg')
            thresh(i) = d_m - (steps(i)*d_std);
        end
    
        %Generate spike times from threshold crossing
        idx = 2:length(data);
        spikeTimes{i} = find(data(idx) > thresh(i) & data(idx-1) < thresh(i))./fs; %up-crossings only
    
        %Populate cell file from the given info
        cell = [];
        n = regexp(keys{i}, '_', 'split');
        cell.birdname = n{1};
        cell.filename = keys{i};
        cell.channel = chan;
        cell.cell_no =  1;
        cell.filenum = str2double(n{2});
        cell.spikes = [];
        cell.noise = [];
        cell.spont = spikeTimes{i};
    
        daystart = 10; %day of month to start on.
        fileHour(i) =  24*(str2double(n{5}) - daystart) + str2double(n{6}) + (str2double(n{7})/60);
        
        %Copy to the exportable structure
        presentCell = [presentCell; cell];
        
        % Parse out the excluded regions of the file based on annotation record
        %Create the segments to exclude (based on annotation record) and a pre/post buffer
        buffer = 150/1000; %segmenting buffer in seconds (250ms)
        cutouts = [(elements{i}.segFileStartTimes-buffer)', (elements{i}.segFileEndTimes+buffer)'];

        %Join overlapping cutouts
        if size(cutouts,1) > 1
            idx = 2:size(cutouts,1);
            o = find((cutouts(idx-1,2) - cutouts(idx,1)) < 0.010); %eliminate gaps of less then 10ms
            t = NaN(size(o,1),2);
            t(:,2) = cutouts(o,2);
            t(1,1) = cutouts(1,1);
            t(2:end,1) = cutouts((o(1:end-1)+1),1);
            cutouts = t;
        else
            a = 1;
        end
        
        %Parse out spikes excluded by the annotaion cutouts
        spikes = spikeTimes{i};
        for j = 1:size(cutouts,1)
            spikes = spikes(spikes<=cutouts(j,1) | spikes>=cutouts(j,2));
        end
        parseSpikes{i} = spikes;
        numSpikes(i) = length(spikes);
        parselength(i) = (length(data)/fs) - sum(abs(diff(cutouts')));
    end

end

spikeRate = numSpikes./parselength;

%Plot data results
h = figure (1); clf
cla; hold on
plot(fileHour, spikeRate, 'ok') %plot spike rate
axis tight; y = ylim;
line([12, 12], [0, y(2)], 'Color', 'r') %Lesion time
line([23.5, 23.5], [0, y(2)], 'Color', 'c') %Lights out
line([34, 34], [0, y(2)], 'Color', 'c') % Lights on
line([46.5, 46.5], [0, y(2)], 'Color', 'c') % Lights out
ylabel(['Spike rate (sp/sec)'])
xlabel(['Time (hours)'])
set(gca,'Box', 'off', 'TickDir', 'out')
hold off
title(['Spontaneous Spiking Rate before and after Nif Lesion'])


%Save the cell file
save([annotName(1:(end-19)) 'spon' num2str(chan) '_cell.mat'], 'presentCell')

%Save the spike rate and timing data
save([annotName(1:(end-19)) 'spon' num2str(chan) '_spikeRate.mat'], 'fileHour', 'spikeRate', 'h')




