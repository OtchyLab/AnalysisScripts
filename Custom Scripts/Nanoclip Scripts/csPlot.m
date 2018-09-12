function csPlot
%Specify the time bounds to plot in this run through
%push = 6;
p1 = [171632, 171706;...
    171708, 171729;...
    171734, 171807;...
    171911, 171920;...
    171930, 171937;...
    171939, 172159];
% sTime = 171436; eTime = 171441;
% sTime = 171526; eTime = 171533;
% sTime = 171548; eTime = 171555;
% sTime = 171632; eTime = 171706;
% sTime = 171708; eTime = 171729;
% sTime = 171734; eTime = 171807;
% sTime = 171911; eTime = 171920;
% sTime = 171930; eTime = 171937;
% sTime = 171939; eTime = 172159;

%List all files in the recordings directory
% homeP = 'C:\Users\Tim\Desktop\LLR32';
homeP = 'V:\SongbirdData\LLR32';
filelist = dir([homeP, filesep, '*.mat']);

%Strip out the time to figureout which to process this time around
numFiles = numel(filelist);


tStamp = [];
for i = 1:numFiles
    s = filelist(i).name;
    sp = regexp(s, '_', 'split');
    sps = sp{end}(1:6);
    tStamp(i) = str2double(sps);
end

tdtPCA = [];
stimID = [];
for push = 1:6
    sTime = p1(push, 1); eTime = p1(push, 2);
    idx = find(tStamp>=sTime & tStamp<=eTime);
    
    %Compile the ID'd files into a single container
    numFiles = numel(idx);
    for i = 1:numFiles
        fullName = [homeP, filesep, filelist(idx(i)).name];
        load(fullName);
        
        %Collapse multiple files into a single data packet
        if i == 1
            tdt = data.tdt.response;
            stimStruct = data.stim.prepulse_us;
            alignTime = data.tdt.times_aligned;
        elseif stimStruct == data.stim.prepulse_us
            tdt = cat(1, tdt, data.tdt.response);
        end
        
    end
    tdtPCA = cat(1, tdtPCA, tdt); %contains all of the stimtrials
    id = push.*ones(size(tdt,1),1);
    stimID = [stimID; id];
    
    %Run the plotting routine
    %plotIt(tdt, alignTime, push)
    
end

%Run the PCA routine
% pcaIt(tdtPCA, alignTime, stimID)

%Run the PCA vs time routine
pcavstimeIt(tdtPCA, alignTime, stimID)

%Run the single-stim/all-chans plots
% allChanIt(tdtPCA, alignTime, stimID)

function pcavstimeIt(tdt, alignTime, stimID)
%Plotting all trials of a current steering trial
trials = size(tdt, 1);
chans = size(tdt, 3);

col = {'b', 'r', 'g', 'k', 'm', 'c'};

t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with

%Zero-align all response traces
resp = [];
xs = [-1, 4];
ys = [-210, 210];
for i = 1:chans
    %Find voltage offset from the pre-stim recording
    delta = mean(tdt(:,1:75,i), 2);
    
    %Zero-align all trials
    resp(:,:,i) = tdt(:,:,i) - delta;
end

%Subtract common mode
for i = 1:chans
    
    % Define averaging mask
    mask = true(1, chans);
    mask(i) = false;
    
    %Process trials separately
    for j = 1:trials
        %Define common mode
        commonMode = mean(resp(j, :, mask), 3);
        
        %Subtract common from signal
        respSub(j,:,i) = resp(j,:,i) - commonMode;
%         respSub(j,:,i) = resp(j,:,i);
        
    end
end
resp = respSub; %copy to main var (simplifying coding)


%Set the response times to analyze
preblock =  [-0.0015, -0.00005];  % 1ms pre
postblock = [0.00075,  0.00575]; % 5ms post
preIdx = alignTime >= preblock(1) & alignTime <= preblock(2);
postIdx = alignTime >= postblock(1) & alignTime <= postblock(2);
idx = preIdx | postIdx;
time = [alignTime(preIdx), NaN, alignTime(postIdx)];
preLength = numel(alignTime(preIdx));

%Concatenate channels across trials
snipSet = [];
for i = 1:trials
    %Extract a single stim trial
    S = squeeze(resp(i,:,:));
    
    %Snip out the Stimulation artifact
    S = S(idx,:)';
    
    %Add to the PCA list
    snipSet = [snipSet, S];

end


%Compute PCA over channel responses
%Normalize (z-score) each of the features over all snips and timepoints
zSnips = zscore(snipSet');
    
%PCA analysis
[coeff, score, latent, ~, explained] = pca(zSnips);

%Unwrap the scores matrix into trials again
snipSize = numel(alignTime(idx));
sc = score';
for i = 1:trials
    starts = snipSize*(i-1)+1;
    ends = snipSize*(i);
    pcaCube(i,:,:) = sc(:, starts:ends);
end

% %Calc mean and sem over traces
stimTypes = unique(stimID);
mPC1 = []; sPC1 = [];
mPC2 = []; sPC2 = [];
mPC3 = []; sPC3 = [];
for i = 1:numel(stimTypes)
    idx = stimID == stimTypes(i);
    
    %Mean and STD of each PC for each stim pattern
    mPC1(i,:) = mean(pcaCube(idx,1,:),1);
    mPC2(i,:) = mean(pcaCube(idx,2,:),1);
    mPC3(i,:) = mean(pcaCube(idx,3,:),1);
    
    sPC1(i,:) = std(pcaCube(idx,1,:),1,1);
    sPC2(i,:) = std(pcaCube(idx,2,:),1,1);
    sPC3(i,:) = std(pcaCube(idx,3,:),1,1);
    
end

%Plot the individual PCA traces by stim type
figure(95); clf

%PC1
subplot(3,1,1)
mPC = [mPC1(:,1:preLength), NaN(6,1), mPC1(:,preLength+1:end)];
sPC = [sPC1(:,1:preLength), NaN(6,1), sPC1(:,preLength+1:end)];
for i = 1:numel(stimTypes)
    %Trim out the artifact period
    plot(time*1000, smooth(mPC(i,:)', 5)', 'Color', colors(i,:), 'LineWidth', 1); hold on
%     shadedErrorBar(time*1000, smooth(mPC(i,:)', 5)', smooth(sPC(i,:)', 5)', col(i), 1); hold on
end
xlim([-1, 4]); ylim([-1, 8]);
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-1:4],'YTick', [-1:3:8])
ylabel('PC1')

%PC2
subplot(3,1,2)
mPC = [mPC2(:,1:preLength), NaN(6,1), mPC2(:,preLength+1:end)];
sPC = [sPC2(:,1:preLength), NaN(6,1), sPC2(:,preLength+1:end)];
for i = 1:numel(stimTypes)
    %Trim out the artifact period
    plot(time*1000, smooth(mPC(i,:)', 5)', 'Color', colors(i,:), 'LineWidth', 1); hold on
%     shadedErrorBar(time*1000, smooth(mPC(i,:)', 5)', smooth(sPC(i,:)', 5)', col(i), 1); hold on
end
xlim([-1, 4]); ylim([-5, 3]);
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-1:4],'YTick', [-5:4:3])
ylabel('PC2')

%PC3
subplot(3,1,3)
mPC = [mPC3(:,1:preLength), NaN(6,1), mPC3(:,preLength+1:end)];
sPC = [sPC3(:,1:preLength), NaN(6,1), sPC3(:,preLength+1:end)];
for i = 1:numel(stimTypes)
    %Trim out the artifact period
    plot(time*1000, smooth(mPC(i,:)', 5)', 'Color', colors(i,:), 'LineWidth', 1); hold on
%     shadedErrorBar(time*1000, smooth(mPC(i,:)', 5)', smooth(sPC(i,:)', 5)', col(i), 1); hold on
end
xlim([-1, 4]); ylim([-0.5, 1]);
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-1:4],'YTick', [-0.5:0.5:1])
ylabel('PC3')
xlabel('Time (ms)')

set(gcf, 'Units', 'Inches', 'Position', [2.75, 6.25, 3.25, 6.5])


function allChanIt(tdt, alignTime, stimID)
%Plotting all trials of a current steering trial
trials = size(tdt, 1);
chans = size(tdt, 3);
stims = numel(unique(stimID));

t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with
ys = [-210, 210];
col = {'b', 'r', 'g', 'k', 'm', 'c'};

%Zero-align all response traces
resp = [];
for i = 1:chans
    %Find voltage offset from the pre-stim recording
    delta = mean(tdt(:,1:75,i), 2);
    
    %Zero-align all trials
    resp(:,:,i) = tdt(:,:,i) - delta;
end

%Subtract common mode
for i = 1:chans
    
    % Define averaging mask
    mask = true(1, chans);
    mask(i) = false;
    
    %Process trials separately
    for j = 1:trials
        %Define common mode
        commonMode = mean(resp(j, :, mask), 3);
        
        %Subtract common from signal
        respSub(j,:,i) = resp(j,:,i) - commonMode;
%         respSub(j,:,i) = resp(j,:,i);
        
    end
end
resp = respSub; %copy to main var (simplifying coding)

%Set the response times to analyze
preblock =  [-0.0015, -0.00005];  % 1ms pre
postblock = [0.00075,  0.00575]; % 5ms post
preIdx = alignTime >= preblock(1) & alignTime <= preblock(2);
postIdx = alignTime >= postblock(1) & alignTime <= postblock(2);
time = [alignTime(preIdx), NaN, alignTime(postIdx)];
%idx = preIdx | postIdx;

%Plot
figure(98); clf
for i = 1:stims
    subplot(3,2,i)
    
    stimMask = stimID == i;
    for j = 1:chans
        %Calc mean and sem over traces
        m = mean(resp(stimMask,:,j)*10e3, 1);
        s = std(resp(stimMask,:,j)*10e3, 1,1)./sqrt(sum(stimMask(:)));
%                 sub = resp(stimMask,:,j)*10e3;
%                 m = squeeze(sub(5,:,:));
        
        %Plot channel common mode subtracted mean and var
%         mTrim = [m(preIdx), NaN, m(postIdx)];
        mTrim = [smooth(m(preIdx)',3)', NaN, smooth(m(postIdx)',3)'];
        sTrim = [s(preIdx), NaN, s(postIdx)];
        plot(time*1000, mTrim,'Color', [0,0,0]+(0.125*j), 'LineWidth', 1.5); hold on
    end
    
    line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
    
    %Format the figure
    axis tight
    xlim([-1,4])
    ylim([-0.225, 0.225])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'YTick', -0.2:0.2:0.2)
    ylabel('Evoked response (mV)')
    xlabel('Time (ms)')
    
end
set(gcf, 'Units', 'Inches', 'Position', [15. 2. 7.5, 11.5])

function pcaIt(tdt, alignTime, stimID)
%Plotting all trials of a current steering trial
trials = size(tdt, 1);
chans = size(tdt, 3);

col = {'b', 'r', 'g', 'k', 'm', 'c'};

t = get(0,'defaultAxesColorOrder');
colors = [t; rand(23,3)]; %30 initial plotting colors to work with

%Zero-align all response traces
resp = [];
xs = [-1, 4];
ys = [-210, 210];
for i = 1:chans
    %Find voltage offset from the pre-stim recording
    delta = mean(tdt(:,1:75,i), 2);
    
    %Zero-align all trials
    resp(:,:,i) = tdt(:,:,i) - delta;
end

%Subtract common mode
for i = 1:chans
    
    % Define averaging mask
    mask = true(1, chans);
    mask(i) = false;
    
    %Process trials separately
    for j = 1:trials
        %Define common mode
        commonMode = mean(resp(j, :, mask), 3);
        
        %Subtract common from signal
        respSub(j,:,i) = resp(j,:,i) - commonMode;
%         respSub(j,:,i) = resp(j,:,i);
        
    end
end
resp = respSub; %copy to main var (simplifying coding)


%Set the response times to analyze
preblock =  [-0.0015, -0.0005];  % 1ms pre
postblock = [0.00075,  0.00575]; % 5ms post
preIdx = alignTime >= preblock(1) & alignTime <= preblock(2);
postIdx = alignTime >= postblock(1) & alignTime <= postblock(2);
idx = preIdx | postIdx;

%Concatenate channels across trials
snipSet = [];
for i = 1:trials
    %Extract a single stim trial
    S = squeeze(resp(i,:,:));
    
    %Snip out the Stimulation artifact
    S = S(idx,:)';
    
    %Add to the PCA list
    snipSet = [snipSet, S];

end


%Compute PCA over channel responses
%Normalize (z-score) each of the features over all snips and timepoints
zSnips = zscore(snipSet');
    
%PCA analysis
[coeff, score, latent, ~, explained] = pca(zSnips);

%Unwrap the scores matrix into trials again
snipSize = numel(alignTime(idx));
sc = score';
for i = 1:trials
    starts = snipSize*(i-1)+1;
    ends = snipSize*(i);
    pcaCube(i,:,:) = sc(:, starts:ends);
end

% %Calc mean and sem over traces
stimTypes = unique(stimID);
mPC1 = []; sPC1 = [];
mPC2 = []; sPC2 = [];
mPC3 = []; sPC3 = [];
for i = 1:numel(stimTypes)
    idx = stimID == stimTypes(i);
    
    %Mean and STD of each PC for each stim pattern
    mPC1(i,:) = mean(pcaCube(idx,1,:),1);
    mPC2(i,:) = mean(pcaCube(idx,2,:),1);
    mPC3(i,:) = mean(pcaCube(idx,3,:),1);
    
    sPC1(i,:) = std(pcaCube(idx,1,:),1,1);
    sPC2(i,:) = std(pcaCube(idx,2,:),1,1);
    sPC3(i,:) = std(pcaCube(idx,3,:),1,1);
    
end


%Plot the individual traces by stim type
figure(105); clf
for i = 1:numel(stimTypes)
    plot(smooth(mPC1(i,:),5), smooth(mPC2(i,:),5),'Color', colors(stimTypes(i),:)); hold on
end

figure(106); clf
for i = 1:numel(stimTypes)
    plot3(smooth(mPC1(i,:),5), smooth(mPC2(i,:),5),smooth(mPC3(i,:),5),'Color', colors(stimTypes(i),:), 'LineWidth', 2); hold on
end
%lims = axis; [az, el] = view;
xlim([-1, 7]); ylim([-4, 2]); zlim([-1, 1]);
view([-43.5,61.2])
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [-1:7],'YTick', [-4:2],'ZTick', [-1:1])
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
grid

function plotIt(tdt, alignTime, push)
%Plotting all trials of a current steering trial
trials = size(tdt, 1);
chans = size(tdt, 3);

col = {'b', 'r', 'g', 'k', 'm', 'c'};

% %Fig w/ subplot per channel -- full response
% figure(100);
% set(gcf, 'Units', 'Inches', 'Position', [1, 1, 4, 9])
%
% %Fig w/ subplot per channel -- common mode subtracted
% figure(101);
% set(gcf, 'Units', 'Inches', 'Position', [6, 1, 4, 9])

%Zero-align all response traces
resp = [];
xs = [-1, 4];
ys = [-210, 210];
for i = 1:chans
    %find offset
    delta = mean(tdt(:,1:75,i), 2);
    
    %Zero-align all trials
    resp(:,:,i) = tdt(:,:,i) - delta;
end

%subtract common mode
for i = 1:chans
    
    % Define averaging mask
    mask = true(1, chans);
    mask(i) = false;
    
    %Process trials separately
    for j = 1:trials
        %Define common mode
        commonMode = mean(resp(j, :, mask), 3);
        
        %Subtract common from signal
        respSub(j,:,i) = resp(j,:,i) - commonMode;
%         respSub(j,:,i) = resp(j,:,i);
        
    end
    
    %Calc mean and sem over traces
    m = mean(respSub(:,:,i)*10e6, 1);
    s = std(respSub(:,:,i)*10e6, 1, 1)./sqrt(trials);
%     s = std(respSub(:,:,i)*10e6, 1, 1);
    
    %Plot channel common mode subtracted mean and var
    figure(102); subplot(3,2,i)
    shadedErrorBar(alignTime*1000, smooth(m, 5), s, col(push), 1); hold on
%     plot(alignTime*1000, smooth(m, 5), col{push}); hold on
    line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
    %Format the figure
    axis tight
    xlim(xs)
    ylim(ys)
    set(gca, 'Box', 'off', 'TickDir', 'out')
    ylabel('Evoked response (\muV)')
    xlabel('Time (ms)')
end


p = 'C:\Users\Tim\Desktop\LLR32';
% savename = [p, filesep, 'Pics', filesep, dataBlock.filename(1:end-3), 'png'];

% pause(1)
% saveas(gcf, savename)




