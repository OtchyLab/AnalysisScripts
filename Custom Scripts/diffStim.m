function diffStim
%Ask user which file to load; check for suffix
[filename, pathname, filterindex] = uigetfile('C:\Users\Tim\Desktop\LLR32\*.mat', 'Pick a files to process', 'MultiSelect', 'on');
numFiles = numel(filename);
dataBlock = [];
for i = 1:numFiles
    fullName = [pathname, filesep, filename{i}];
%     fullName = [pathname, filesep, filename];
    load(fullName);
    dataBlock = data;
    
    plotIt(dataBlock)
end


function plotIt(dataBlock)
%Plotting all trials of a current steering trial


%Plotting sepoarately all trials of a single sTDT recording
trials = size(dataBlock.tdt.response, 1);
chans = size(dataBlock.tdt.response, 3);

col = {'b', 'r', 'g', 'k', 'm', 'c'};

%Single figure with all of the channels in subplots
% figure(6969); clf
% set(gcf, 'Units', 'Inches', 'Position', [0.5, 0.5, 3.5, 9])
% 
% figure(6868); clf
% set(gcf, 'Units', 'Inches', 'Position', [4.5, 0.5, 3.5, 9])

f = figure(6767); clf
set(gcf, 'Units', 'Inches', 'Position', [8, 6, 8.25, 2.75])

%Zero-align all response traces
resp = [];
xs = [-2, 4];
ys = [-150, 150];
for i = 1:chans
    %find offset
    delta = mean(dataBlock.tdt.response(:,1:75,i), 2);
    
    %Zero-align all trials
    resp(:,:,i) = dataBlock.tdt.response(:,:,i) - delta;
end

%subtract common mode
for i = 1:chans
    %Select figure/subplot
    %figure(6969); subplot(chans, 1, i)
    
    % Define averaging mask
    mask = true(1, chans);
    mask(i) = false;
    
    %Process trials separately
    for j = 1:trials
        %Define common mode
        commonMode = mean(resp(j, :, mask), 3);
        
        %Subtract common from signal
        respSub(j,:,i) = resp(j,:,i) - commonMode;
        
        %Plot
        %plot(dataBlock.tdt.times_aligned*1000, respSub(j,:,i)*10e6); hold on
    end
    %line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
    
%     %Format the figure
%     if i ~= chans
%         set(gca, 'XTickLabels', '');
%     else
%         xlabel(['Time (ms)'])
%     end
%     axis tight
%     xlim(xs)
%     ylim(ys)
%     set(gca, 'Box', 'off', 'TickDir', 'out')
%     ylabel(['Channel ' num2str(i), ' (\muV)'])
    
    %Select figure/subplot
    %figure(6868); subplot(chans, 1, i)
    
    m = mean(respSub(:,:,i)*10e6, 1);
%     s = std(respSub(:,:,i)*10e3, 1, 1);
    s = std(respSub(:,:,i)*10e6, 1, 1)./sqrt(trials);
    %shadedErrorBar(dataBlock.tdt.times_aligned*1000, m, s, 'k', 1)
    %line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');
    
%     %Format the figure
%     if i ~= chans
%         set(gca, 'XTickLabels', '');
%     else
%         xlabel(['Time (ms)'])
%     end
%     axis tight
%     xlim(xs)
%     ylim(ys)
%     set(gca, 'Box', 'off', 'TickDir', 'out')
%     ylabel(['Channel ' num2str(i), ' (\muV)'])
    
    %Select figure
%     figure(6767);
    shadedErrorBar(dataBlock.tdt.times_aligned*1000, m, s, col(i), 1); hold on
    
    %Format the figure
    axis tight
    xlim(xs)
    ylim(ys)
    set(gca, 'Box', 'off', 'TickDir', 'out')
    ylabel('Evoked response (\muV)')
    xlabel('Time (ms)')
end
line([0, 0], ys, 'Color', 'r', 'LineStyle', ':');

p = 'C:\Users\Tim\Desktop\LLR32';
savename = [p, filesep, 'Pics', filesep, dataBlock.filename(1:end-3), 'jpeg'];
%savefig(f, savename)
saveas(f, savename, 'jpeg') 

close(f)








