function acute_sweep_analysis
%This function will take as input a user-specified source folder and range
%of files and do some plotting and summary stats on the selected files.

clear

%Request input from the user on which folder to process
% folder = '/Users/tim/Desktop/Acute Recordings/LLR20';
folder = '/Users/tim/Desktop/Acute Recordings/LLR32';
% folder = '/Users/tim/Desktop/Acute Recordings/Rd40';

%A few current sweeps from LLR20
% start_time = 145924;
% end_time = 150220;

start_time = 154457;
end_time = 154628;

%Processing constants
stimMax = 95;
stimMin = 5;
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
    
    %Cycle through the channels
    for j = 1:6
        %Set plotting location
        subplot(2, 3, j)
        caxis([stimMin, stimMax])
        
        %Subtract the mean to center the resposne
        delta = (mean(dataBlock.tdt.response(:,1:75,j),2));
        mean_sub_responses = (dataBlock.tdt.response(1:1:end,:,j)-delta(1:1:end))*10e3; 
        
        %Handle each trial 
        for k = 1:num_trials
            %Screen for noise on the pre-stim signal
            rec = mean_sub_responses(k,:);
            snip = rec(tMask);
            resp = rec(rMask);
            if min(snip) > -0.25 && max(snip) < 0.25
%              if true
                %Plot it
                plot(times, rec, 'Color', col); hold on 
                
                %Calculate Vpp
                vpp = max(resp) - min(resp);
                
                %Shuffle info into carry-over function
                summary = [summary; j, stim, vpp];
            end
        end
        
        %Format the figure
        axis tight
        xlim([-2, 6])
        ylim([-2, 3])
        set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', -2:1:2)
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
    
    xlim([stimMin, stimMax]); ylim([0, 5])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', 0:2.5:5)
    xlabel('Current (uA)'); ylabel('Vpp (mV)')

end

%Finished cooking
disp('...')
disp('Finished processing dataset')
disp('...')





