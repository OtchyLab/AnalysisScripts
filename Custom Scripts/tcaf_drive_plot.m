%reads the tCAF log file (NOT the detect file!) and makes a pretty heatmap
%of the target durations

%% Reset what you will...
clear

%Set the file to load and process
% mother = '/Users/Tim/Desktop';
mother = 'C:\Users\Tim\Desktop';
file = 'LW39_Drive1';
loc = [mother filesep file];

%Load the drive file and parse ascii data to a cell array
fileID = fopen(loc);
C = textscan(fileID,'%s %s %f %u8 %f %u8 %f %u8 %f %f %f %f %f %s');
fclose(fileID);
%MM/dd/yyyy HH:mm:ss.SSS
%Copy out the data I care about to meaningful variables
dates = C{1};
times = C{2};
ltimes = datetime([char(dates), char(times)],'InputFormat','MM/dd/yyyyHH:mm:ss.SSS');
dur = C{3};
jam = C{4};
ave = C{7};
catches = C{8};
wavnames = C{14};

%% Plot heatmap of durations
% Filter parameters
window = 10;
factor = 2;

%Other parameters
maxInterval = 200;
trialBinSize = 5;
intertapBinSize = 2;
upperScale = 0.125;

%Assemble the heatmap
TargetDur = dur;
intertapMatrix = zeros(length(TargetDur), maxInterval);
convMatrix = ones(trialBinSize, intertapBinSize);
for m = 1:length(TargetDur)
    if (TargetDur(m) > 0) && (TargetDur(m) < maxInterval)
        intertapMatrix(m,round(TargetDur(m))) = 1;
    end
end
matrix = conv2(intertapMatrix, convMatrix, 'same');

% Filter
hfilt = ones(window*factor, window)/(factor*window^2);
matrix = imfilter(matrix, hfilt, 'replicate');
    
% Normalize in each rendition-bin
for m = 1:size(matrix,1)
    matrix(m,:) = matrix(m,:)/sum(matrix(m,:));
end

%Find the day starts
ds = datenum(char(dates), 'MM/dd/yyyy');
idx = find(diff(ds)==1);
idx = [1; idx];

%Plot to figure
h = figure(1); clf

imagesc(matrix', [0 upperScale]); colormap(jet)
hold on
axis xy

for i = idx
    l = plot([i, i], [0, maxInterval], '-k', 'LineWidth', 1.5);
end

set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Rendition', 'FontSize', 18)
ylabel('Target Duration (ms)', 'FontSize', 18)
set(gcf, 'Units', 'Inches', 'Position', [5, 4.5, 9, 5.5], 'LineWidth', 3, 'FontSize', 18)


