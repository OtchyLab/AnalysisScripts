function csPlot
%Specify the time bounds to plot in this run through
%push = 6;
p1 = [171632, 171706;...
    171708, 171729;...
    171734, 171807;...
    171911, 171920;...
    171930, 171937;...
    171939, 172159];
% % sTime = 171436; eTime = 171441;
% % sTime = 171526; eTime = 171533;
% % sTime = 171548; eTime = 171555;
% sTime = 171632; eTime = 171706;
% sTime = 171708; eTime = 171729;
% sTime = 171734; eTime = 171807;
% sTime = 171911; eTime = 171920;
% sTime = 171930; eTime = 171937;
% sTime = 171939; eTime = 172159;

%List all files in the recordings directory
homeP = 'C:\Users\Tim\Desktop\LLR32';
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
    
    %Run the plotting routine
    plotIt(tdt, alignTime, push)
end

function plotIt (tdt, alignTime, push)
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
        
    end
    
    %     %Calc mean and sem over traces
    %     m = mean(resp(:,:,i)*10e6, 1);
    %     s = std(resp(:,:,i)*10e6, 1, 1)./sqrt(trials);
    %
    %     %Plot channel common mode subtracted mean and var
    %     figure(i); %subplot(chans,1,i)
    %     shadedErrorBar(dataBlock.tdt.times_aligned*1000, m, s, col(i), 1); hold on
    %     plot(alignTime*1000, smooth(m, 5), col{i}); hold on
    
    %     %Format the figure
    %     axis tight
    %     xlim(xs)
    %     ylim(ys)
    %     set(gca, 'Box', 'off', 'TickDir', 'out')
    %     ylabel('Evoked response (\muV)')
    %     xlabel('Time (ms)')
    
    %Calc mean and sem over traces
    m = mean(respSub(:,:,i)*10e6, 1);
    s = std(respSub(:,:,i)*10e6, 1, 1)./sqrt(trials);
%     s = std(respSub(:,:,i)*10e6, 1, 1);
    
    %Plot channel common mode subtracted mean and var
    figure(101); subplot(3,2,i)
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




