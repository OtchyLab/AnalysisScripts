function acute_sweep_analysis_thinCMS
%This function will take as input a user-specified source folder and range
%of files and do some plotting and summary stats on the selected files.
%
%
% Created by TMO 9/30/20

%Data selection

%A few current sweeps from LLR20
folder = '/Users/tim/Desktop/Acute Recordings/LLR20';
start_time = 145924;
end_time = 150220;

% start_time = 145914;
% end_time = 145937;

%A few current sweeps from Rd40
% folder = '/Users/tim/Desktop/Acute Recordings/Rd40';
% start_time = 162404;
% end_time = 163425;

%A few current sweeps from LLR32
% folder = '/Users/tim/Desktop/Acute Recordings/LLR32';
% start_time = 153837;
% end_time = 154019;

%Processing constants/params
stimRng = [0, 110]; %uA
respRng = [0.75, 6]; %ms
select = [1, 2, 5, 6]; %chans
select = 1:6;

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

%Extract trial data from selected files
[trial_stims, trial_times, trial_responses] = extract_trials(folder, sub_files); 

%Compute the metrics we care about
trial_vpps = calc_response_metrics(trial_responses, trial_times, respRng);

%Plot binned sweeps for each selected channel
h1 = plot_binned_sweeps(trial_responses, trial_times, trial_stims, select, stimRng);

%Plot the in-out curve for each channel
h2 = plot_in_out(trial_stims, trial_vpps, select, stimRng);

%Plot the binned in-out curve for each channel
h3 = plot_binned_in_out(trial_stims, trial_vpps, select, stimRng);

%Finished cooking
disp('...')
disp(['Finished processing dataset'])
disp('...')

function [trial_stims, trial_times, trial_responses] = extract_trials(folder, sub_files)
%Sequentially process selected files
trial_stims = [];
trial_times = [];
trial_responses = [];

%  for i = 1:15
for i = 1:numel(sub_files)
    %Load the current file
    cur_file = sub_files(i).name;
    load([folder, filesep, cur_file]);
    dataBlock = data;

    %Extract params from the current file
    stim = dataBlock.stim.current_uA;
    num_trials = size(dataBlock.tdt.response, 1);
    times = dataBlock.tdt.times_aligned*1000; 
    
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
        % Define averaging mask
        mask = true(1, 6);
        mask(j) = false;
        
        %Process trials separately
        for k = 1:num_trials
            %Define common mode
            commonMode = mean(demean_responses(k, :, mask), 3);
            
            %Subtract common from signal
            CMS_responses(k,:,j) = demean_responses(k,:,j) - commonMode;
            
        end
    end
    
  
    %Copy out data to the carrying variables
    
    many_stim = repmat(stim,num_trials, 1);
    trial_stims = [trial_stims; many_stim];
    
    many_times = repmat(times,num_trials, 1);
    trial_times = [trial_times; many_times];
    
    trial_responses = [trial_responses; CMS_responses];
end    

function trial_vpps = calc_response_metrics(responses, times, respRng)

%Define output
trial_vpps = [];

% Mask for later analysis
rMask = times(1,:) >= respRng(1) & times(1,:) < respRng(2);

%Cycle through the trials and channels to extract the quantities
num_trials = size(responses, 1);
num_chans = size(responses, 3);
for i = 1:num_trials
   for j = 1:num_chans
        %Select an individual recording
        rec = responses(i,:,j);
        roi = rec(rMask);
        
        %Calc and store VPP
        trial_vpps(i,j) = max(roi) - min(roi);
   end
end

function h = plot_in_out(stims, vpps, select, stimRng)
%Plot the in-out cloud for each channel
h = figure(2); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 10, 8])


for i = 1:numel(select)
        %Select active subplot
        subplot(2, 3, i)
        
        %Set channel
        chan = select(i);
    
        %Scatter plot it all
        scatter(stims, vpps(:, chan), '.k')
    
        xlim([stimRng(1)-10, stimRng(2)+10]); ylim([0, 1])
        set(gca, 'Box', 'off', 'TickDir', 'out')%, 'XTick', 0:50:stimRng(1))
        xlabel('Current (uA)'); ylabel('Vpp (mV)')

end

function h = plot_binned_in_out(stims, vpps, select, stimRng)
%Plot the mean in-out curve for channels
h = figure(3); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 10, 8])

%Binning stuff
step = 10;
edges = stimRng(1):step:stimRng(2);
center = edges - (step/2);
bins = discretize(stims,edges);
unique_bins = sort(unique(bins));

%Cycle through the channels
vpps_mean = [];
vpps_std = [];
for i = 1:numel(select)
    chan = select(i);
    
    %Step through the bins
    for j = 1:numel(unique_bins) 
        
        mask = bins == unique_bins(j);
        %Calc mean and std for all channel trials in that bin
        vpps_mean(i, j) = mean(vpps(mask, chan), 1);
        vpps_std(i, j) = std(vpps(mask, chan), 1);
    end
    
    %Plot this channel
    subplot(2, 3, i)
    errorbar(center(unique_bins), vpps_mean(i, :), vpps_std(i, :), '.-')
    
    xlim([stimRng(1)-10, stimRng(2)+10]); ylim([0, 1])
    set(gca, 'Box', 'off', 'TickDir', 'out')%, 'XTick', 0:50:stimRng(1))
    xlabel('Current (uA)'); ylabel('Vpp (mV)')
    
end

function h = plot_binned_sweeps(responses, times, stims, select, stimRng)

%Setup the figure
h = figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 10, 8])
clr = parula(100);

%Binning stuff
step = 10;
edges = stimRng(1):step:stimRng(2);
center = edges - (step/2);
bins = discretize(stims,edges);
unique_bins = sort(unique(bins));

%Step through the bins
for i = 1:numel(unique_bins)
    cur_bin = unique_bins(i);
    cur_stim = center(cur_bin);
    
    %Figure out trace color based on the stim current
    ramp = linspace(stimRng(1), stimRng(2), 100);
    a = interp1(ramp, clr(:,1), max([0, cur_stim]));
    b = interp1(ramp, clr(:,2), max([0, cur_stim]));
    c = interp1(ramp, clr(:,3), max([0, cur_stim]));
    col = [a, b, c];
    
    %Cycle through the channels
    mask = bins == cur_bin;
    
    for j = 1:numel(select)
        subplot(2, 3, j)
        chan = select(j);
        
        %Calc mean and std for all channel trials in that bin
        resp_mean = mean(responses(mask,:,chan), 1);
        resp_std = std(responses(mask,:,chan), 1);
        
        %Plot it
        t = times(1,:);
        t(t>0.02 & t<0.75) = NaN;
        plot(t, resp_mean, 'Color', col, 'LineWidth', 1.5); hold on
        
                
        %Format the figure
        axis tight
        xlim([-2, 6])
        ylim([-0.6, 0.6])
        set(gca, 'Box', 'off', 'TickDir', 'out')%, 'YTick', -0.5:0.5:0.5)
        xlabel('Time (ms)'); ylabel(['Evoked Response (mV)'])
    end
end


 
     

  
  



