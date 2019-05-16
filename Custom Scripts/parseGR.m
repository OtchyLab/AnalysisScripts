%Script for parsing avis from the 1P2C scope
%
%Written by TMO 05/14/2019

%Load raw interleaved AVI from file w/ MATLAB native video loader
video = []; ind = 0;
v = VideoReader('C:\Users\Tim\Downloads\04.26.2019Interleaved.avi');
while hasFrame(v)
    ind = ind +1;
    video(ind,:,:) = readFrame(v);
end

%Correlation frame checking
template = squeeze(video(1,:,:));
template = template(:);
xc = 1;
for i = 2:size(video,1)
    %Vectorize frame
   curFrame =  squeeze(video(i,:,:));
   curFrame = curFrame(:);
    
   %Xcorr
   [xc(i), ~] = corr(template, curFrame);
end

%Fit 2-component gaussian model to the correlation data
gm = fitgmdist(xc',2);
thresh = mean(gm.mu); %Could do some error checking here...

%Get frame indices based on frame correlation
gInd = xc>thresh;
rInd = xc<thresh;

%Plot frame correlation and thresholding 
figure(667); clf

plot(xc); hold on
line([1, numel(xc)], [thresh, thresh], 'Color', 'r');
xlabel('Frame #'); ylabel('\rho')
set(gca, 'Box', 'off', 'TickDir', 'out')
set(gcf, 'Units', 'Inches', 'Position', [10, 10, 11.5, 3.25])

%Parse interleaved frames into channels
gFrames = video(gInd,:,:);
rFrames = video(rInd,:,:);

%Average over frames to reduce sensor/thermal noise
gMean = squeeze(mean(gFrames,1));
rMean = squeeze(mean(rFrames,1));

%Plot mean images
figure(666); clf

subplot(1,3,1)
imagesc(gMean, [0, 255])
colormap('gray')
title('Mean Green Frame')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

subplot(1,3,2)
imagesc(rMean, [0, 255])
title('Mean Red Frame')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

subplot(1,3,3)
imagesc(gMean-rMean, [0, 255])
title('Mean Frame Delta')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [5.25, 9.5, 14, 3.25])


%Psuedo colored image plot
%Scaling?
rScale = 255/max(rMean(:));
gScale = 255/max(gMean(:));

% Create an all black frame.
allBlack = zeros(size(rMean), 'uint8');

% Create RGB channels.
redChannel = rMean.*rScale; % Red channel
greenChannel = gMean.*gScale; % Green channel
blueChannel = allBlack; % Blue channel

% Create color versions of the two color channels.
just_red = cat(3, redChannel, allBlack, allBlack);
just_green = cat(3, allBlack, greenChannel, allBlack);

% Recombine the individual color channels to create the original RGB image again.
recombinedRGBImage = cat(3, redChannel, greenChannel, blueChannel);

%Plot mean images
figure(668); clf

subplot(1,3,1)
imagesc(just_green, [0, 255])
colormap('gray')
title('Mean Green Frame')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

subplot(1,3,2)
imagesc(just_red, [0, 255])
title('Mean Red Frame')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

subplot(1,3,3)
imagesc(recombinedRGBImage, [0, 255])
title('Mean Frame Merge')
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])
set(gcf, 'Units', 'Inches', 'Position', [5.25, 9.5, 14, 3.25])

