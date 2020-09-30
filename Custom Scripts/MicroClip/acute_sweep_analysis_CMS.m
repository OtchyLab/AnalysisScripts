%%function acute_sweep_analysis
%This function will take as input a user-specified source folder and range
%of files and do some plotting and summary stats on the selected files.

clear

%Request input from the user on which folder to process
folder = '/Users/tim/Desktop/Acute Recordings/LLR20';
% folder = '/Users/tim/Desktop/Acute Recordings/Rd40';
% folder = '/Users/tim/Desktop/Acute Recordings/LLR32';

%A few current sweeps from LLR20
start_time = 145924;
end_time = 150220;

% start_time = 145914;
% end_time = 145937;

%A few current sweeps from Rd40
% start_time = 162404;
% end_time = 163425;

%A few current sweeps from LLR32
% start_time = 153837;
% end_time = 154019;

%Processing constants
stimMax = 110;
stimMin = 10;
respMin = 0.75;
respMax = 6;
summary = [];

%Set figure
figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 16.25, 8])
clr = parula(100);

%Get the folder file contents and check for empties
files = dir([folder, filesep, 'stim*.mat']);
if isempty(files)
    disp('Nothing in the folder, bro... try again.')
    return
end

%Extract timestamp from the filenames
tStamp = [];
for i = 1:numel(files)
    s = files(i).name;
    sp = regexp(s, '_', 'split');
    sps = sp{end}(1:6);
    tStamp(i) = str2double(sps);
end

%Derive selected files from masking on time bounds
mask = (tStamp >= start_time) & (tStamp <= end_time);
sub_files = files(mask);

%Sequentially process selected files

% for i = 1:numel(sub_files)
for i = 1:numel(sub_files)
    %Load the current file
    cur_file = sub_files(i).name;
    load([folder, filesep, cur_file]);
    dataBlock = data;

    %Extract the stimulation current and map to plotting color
    stim = dataBlock.stim.current_uA;
    a = interp1(linspace(stimMin,stimMax,100), clr(:,1), stim);
    b = interp1(linspace(stimMin,stimMax,100), clr(:,2), stim);
    c = interp1(linspace(stimMin,stimMax,100), clr(:,3), stim);
    col = [a, b, c];
        
    %Find number of trials in this file
    num_trials = size(dataBlock.tdt.response, 1);
    times = dataBlock.tdt.times_aligned*1000; 
    tMask = (times >= -2 & times < -0.05) | (times >= 5 & times < 25);
    rMask = times >= respMin & times < respMax;
    
    %Cycle through the channels to demean
    demean_responses = [];
    CMS_responses = [];
    for j = 1:6
        %Subtract the mean to center the resposne
        delta = (mean(dataBlock.tdt.response(:,1:75,j),2));
        demean_responses(:,:,j) = (dataBlock.tdt.response(:,:,j)-delta)*10e3; 
    end
    
    %Cycle through channels to calculate the CMS
    for j = 1:6
        %Set plotting location
        subplot(2, 3, j)
        caxis([stimMin, stimMax])
        
        % Define averaging mask
        mask = true(1, 6);
        mask(j) = false;
        
        %Process trials separately
        for k = 1:num_trials
            %Define common mode
            commonMode = mean(demean_responses(k, :, mask), 3);
            
            %Subtract common from signal
            CMS_responses(k,:,j) = demean_responses(k,:,j) - commonMode;

            %Plot it
            rec = CMS_responses(k,:,j);
            resp = rec(rMask);
            plot(times, rec, 'Color', col); hold on
            
            
            %Calculate Vpp
            vpp = max(resp) - min(resp);
            
            %Shuffle info into carry-over function
            summary = [summary; j, stim, vpp];
 
        end
        
        %Format the figure
        axis tight
        xlim([-2, 6])
        %ylim([-0.6, 0.6])
        set(gca, 'Box', 'off', 'TickDir', 'out')%, 'YTick', -0.5:0.5:0.5)
        xlabel('Time (ms)'); ylabel(['Evoked Response (mV)'])
        
    end
    
    
end

%Plot the in-out curve for each channel
figure(2); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 16.25, 8])
for i = 1:6
    %Set plotting location
    subplot(2, 3, i)
    
    %Select everthing for the current channel
    chan_mask = summary(:,1) == i;
    sub_summary = summary(chan_mask,:);
    
    scatter(sub_summary(:,2), sub_summary(:,3), '.k')
    
    xlim([stimMin, stimMax]); %ylim([0, 1])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.5:1)
    xlabel('Current (uA)'); ylabel('Vpp (mV)')

end

%%
%Plot the mean in-out curver for channels
figure(3); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 16.25, 8])
for i = 1:6
    %Set plotting location
    subplot(2, 3, i)
    
    %Select everthing for the current channel
    chan_mask = summary(:,1) == i;
    sub_summary = summary(chan_mask,:);
    
    stims = sub_summary(:,2);
    vpps = sub_summary(:,3);
    step = 5;
    edges = 10:step:110;
%     edges = 0:step:75;
    center = edges - (step/2);
    bins = discretize(stims,edges);
    unique_bins = sort(unique(bins));
    m_vpp = []; s_vpp = [];
    for j = 1:numel(unique_bins)
       mask_vpp = bins == unique_bins(j);
       m_vpp(j) = mean(vpps(mask_vpp));
       s_vpp(j) = std(vpps(mask_vpp));
    end
    
    errorbar(center(unique_bins), m_vpp, s_vpp, '.-')
    
    xlim([stimMin, stimMax]); %ylim([0, 1])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:0.5:1)
    xlabel('Current (uA)'); ylabel('Vpp (mV)')

end

%Finished cooking
disp('...')
disp(['Finished processing dataset'])
disp('...')





