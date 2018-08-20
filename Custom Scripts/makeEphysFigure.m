%This script is to make an example figure of the song amplitude,
%spectrogram, raw ephys recording, and common mode subtracted recording

%Source file
sourceDir = 'V:\SongbirdData\LR81RY177\2017-09-04';
filename = 'LR81RY177_3587373336_2017_09_04_08_35_36.412.dat';

%Constants
timeS = 2.6; %in seconds from start
timeE = 4.4; %in seconds from start
chan = 1; %ephys channel of interest

%Load the data from file
[channels,fs] = getChannels([sourceDir, filesep, filename]);

%Convert from time to the sample count
startT = timeS * fs;
endT = timeE * fs;

%Trim data to the desired section
channelsTrim = channels(:, startT:endT);

%Audio data filtering
HP_fNorm = 300/(fs/2);
LP_fNorm = 10000/(fs/2);
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

%Apply amplification and filters
filtAudio = filtfilt(BP_b, BP_a, channelsTrim(1,:)-mean(channelsTrim(1,:)));

%Ephys data filtering
HP_fNorm = 300/(fs/2);
LP_fNorm = 15000/(fs/2);
[BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

%Apply amplification and filters
filtEphys = filtfilt(BP_b, BP_a, channelsTrim(chan+1,:)-mean(channelsTrim(chan+1,:)));

%Calculate common mode
CM = mean(channelsTrim(3:5,:),1);



%Plot it all...
figure(99); clf
xs = (1:numel(filtAudio))*1000/fs;

%Amplitude
subplot(4,1,1)
plot(xs, filtAudio, 'b')

axis tight
ylim([-1.2, 1.2]); ylabel('Sound Amplitude')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:500:2000, 'XTickLabels', [])

%Spectrogram
subplot(4,1,2)
displaySpecgramQuick(filtAudio, fs,[0,8500],[],0);

axis tight
xlabel([]); ylabel('Frequency')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTickLabels', [], 'YTickLabels', {'0', '', '', '', '8k'})

%Raw Ephys
subplot(4,1,3)
plot(xs, filtEphys, 'k')

axis tight
xlabel([]); 
ylim([-0.5, 0.5]); ylabel('Raw Activity (V)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:500:2000, 'XTickLabels', [])

%CMS Ephys
subplot(4,1,4)
plot(xs, channelsTrim(chan+1,:) - CM, 'k')

axis tight
xlabel('Time (ms)'); 
ylim([-0.2, 0.2]); ylabel('CMS Activity (V)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 0:500:2000)

%Figure format
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 8])
