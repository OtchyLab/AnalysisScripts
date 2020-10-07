function acute_lido
%This function will take as input a user-specified source folder and range
%of files and do some plotting and summary stats on the selected files.
%
%
% Updated by TMO 9/29/20
clear

%Request input from the user on which folder to process
folder = '/Users/tim/Desktop/Acute Recordings/LLR20';

%A few current sweeps from LLR20
saline_start = 150448;
saline_end = 150505;

lido_start = 150648;
lido_end = 150652;

washout_start = 151746;
washout_end = 151813;

%Processing constants
respMin = 0.75;
respMax = 6;
select = [1, 2, 5, 6];

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
saline_mask = (tStamp >= saline_start) & (tStamp <= saline_end);
lido_mask = (tStamp >= lido_start) & (tStamp <= lido_end);
washout_mask = (tStamp >= washout_start) & (tStamp <= washout_end);

%Data collection
[saline_responses, times] = collect_data(folder, files(saline_mask));
lido_responses = collect_data(folder, files(lido_mask));
washout_responses = collect_data(folder, files(washout_mask));

%Process data
saline_CMS = process_data(saline_responses);
lido_CMS = process_data(lido_responses);
washout_CMS = process_data(washout_responses);

%Plot the data
f1 = figure(1); clf
set(gcf, 'Units', 'Inches', 'Position', [10, 7.75, 10, 8])
plot_data(f1, times, saline_CMS, select, 'k')
plot_data(f1, times, lido_CMS, select, 'r')
plot_data(f1, times, washout_CMS, select, 'g')

%Calculate Vpp
mask = times >= respMin & times <= respMax;
saline_vpp = vpp_data(saline_CMS, mask);
lido_vpp = vpp_data(lido_CMS, mask);
washout_vpp = vpp_data(washout_CMS, mask);



%%%%%%%%%%%%%
%Subfunctions

function vpp = vpp_data(responses, mask)
vpp = [];
for i = 1:6
    mean_trace = mean(responses(:,:,i))*10e3;
    snip = mean_trace(mask);
    vpp = [vpp; max(snip) - min(snip)];
end


function CMS_responses = process_data(responses)

%Demean
for i = 1:6
    %find offset
    delta = (mean(responses(:,1:75,i),2));

    responses(:,:,i) = responses(:,:,i) - delta;
end

%CMS
for i = 1:6
    % Define averaging mask
    mask = true(1, 6);
    mask(i) = false;

    %Process trials separately
    for j = 1:size(responses,1)
        %Define common mode
        commonMode = mean(responses(j, :, mask), 3);
        
        %Subtract common from signal
        CMS_responses(j,:,i) = responses(j,:,i) - commonMode;
    end
end


function plot_data(f, times, responses, select, col)
figure(f)

%Blank artifact
T = times;
preT = find(times < -0.05);
postT = find(times > 0.75);

for j = 1:numel(select)
    i = select(j);
    times = T;
    %Plot the mean trace
    subplot(2, 2, j)
    shadedErrorBar(times(preT), mean(responses(:,preT,i))*10e3, std((responses(:,preT,i))*10e3, 1, 1), col, 1); hold on
    if strcmp(col, 'g')
        times = times-0.2;
        postT = find(times > 0.75);
    end
    shadedErrorBar(times(postT), mean(responses(:,postT,i))*10e3, std((responses(:,postT,i))*10e3, 1, 1), col, 1); hold on
    
    %Format the figure
    axis tight
    xlim([-2, 6])
    ylim([-0.35, 0.35])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', -2:2:6, 'YTick', -0.35:0.35:0.35)
    ylabel('Evoked Response (mV)')
    xlabel('Time (ms)')
    
end


function [responses, times] = collect_data(folder, files)
responses = [];
for i = 1:numel(files)
    cur_file = files(i).name;
    load([folder, filesep, cur_file]);
    
    responses = cat(1, responses, data.tdt.response);
end
times = data.tdt.times_aligned*1000;




















