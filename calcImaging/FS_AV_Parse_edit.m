% this should probably also rename/move the csv files 

function FS_AV_Parse(DIR,varargin)
% FS_AV_Parse

% Parse Data from FreedomScopes
%   Created: 2015/08/02
%   By: WALIII
%   Updated: 2017/09/23 % added gain offset.
%   By: WALIII
%   Updated: 2018/07/09 % modify file saving for use with the talon pipeline
%   By: sadieiq694

% FS_AV_parse will do several things:
%
%   1. Seperate Audio and Video, and place them into .mat files in a .mat
%      directory
%   2. Make spectrogram .gif files for extracted audio, for perusing
%      manually.
%   3. Uses the dependancy extractmedia(), an amazing script created by Nathan Perkins

% Run in the Directory of the .mov files. FS_AVparse should be run first:
% FS_AVparse-->FS_TemplateMatch-->FS_FS_Plot_ROI--> BatchDff2----> FS_Image_ROI--> FS_FS_Plot_ROI

% % navigate to the directory containing the .mov files
cd '/Users/sadiela/Desktop/LW9'
% pull out list of all csv files in the directory
% mov_files = dir('*.mov');
% set target directory
targetDirectory = '/Users/sadiela/Desktop/calc_data';

birdName = getBirdName(pwd);

targetDirectory=strcat(targetDirectory, '/', birdName);
%gif_dir=[pwd,'/gif'];
%error_dir =[pwd,'/error'];

if exist(targetDirectory,'dir'), rmdir(targetDirectory,'s'); end
%if exist(gif_dir,'dir') rmdir(gif_dir,'s'); end
% if exist(error_dir,'dir') rmdir(error_dir,'s'); end

mkdir(targetDirectory); % big directory for the bird data; need another one
% for individual dates
%mkdir(gif_dir);
%mkdir(error_dir);



outlier_flag=0;
if nargin<1 | isempty(DIR), DIR=pwd; end
% you can pass in the directory you want to sort/process

mov_listing=dir(fullfile(DIR,'*.mov'));
mov_listing={mov_listing(:).name};

filenames=mov_listing; % just names of all of the .mov files

csvfiles = dir(fullfile(DIR, '*.csv'));
csvfilenames ={csvfiles(:).name};



disp('Parsing Audio and Video files');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)
    
    % create specific directory based on date of .mov file
    date = filenames{i}(1:10);
    specTarget = strcat(targetDirectory, "/", date);
    if ~exist(specTarget, 'dir'), mkdir(specTarget); end

    [path,file,ext]=fileparts(filenames{i});


	fprintf(1,formatstring,round((i/length(mov_listing))*100));
    FILE = fullfile(DIR,mov_listing{i})
    

f = dir(FILE)
if f.bytes/1000000< 900.000

%Extract Audio and video data...
try
[a_ts, a, v_ts, v] = extractmedia(FILE);
catch
   v1 = VideoReader(FILE);
   k = 1;
    while hasFrame(v1)
    v{k} = readFrame(v1);
    k = k+1;
    end
    v = v';
    a = 0;
    a_ts = 0;
    v_ts = 0;
end
    clear k V1;

else
    
  disp('moving file- to large for batch processing... use FS_AV_Parse(pwd,large)');
  %LargeDir = strcat(path,'/','LargeFiles');
  %movefile(FILE, error_dir)
  continue;
end





% Format VIDEO DATA
[video.width, video.height, video.channels] = size(v{1});
video.times = v_ts% 0.1703*(day-1);

video.nrFramesTotal = size(v,1);
video.FrameRate = 1/mean(diff(v_ts));
for ii = 1: size(v,1)
    est(:,:,:,ii) = v{ii};
end

%
disp('Performing Gain correction')
noise = squeeze(est(:,:,3,:)); % blue channel
noise = (((squeeze(mean(mean(noise(:,1:80,:),1))))));%+(squeeze(mean(mean(noise(1:10,:,:),2)))))/2);

sig = squeeze(est(:,:,2,:)); % green channel
sig = (squeeze(mean(mean(sig(:,1:80,:),1))));
video.gain = noise;
for ii = 1: size(v,1)
     %video.frames(:,:,:,ii) = v{ii}-(noise(ii,:)-min(noise(30:end,:)));
      video.frames(:,:,:,ii) = v{ii}; %-(noise(ii,:)-min(noise(30:end)));
end
clear noise;
clear est;


% Format AUDIO DATA
audio.nrChannels = 1;
audio.bits = 16;
audio.nrFrames = length(a);
audio.data = double(a);
audio.rate = 48000;
audio.TotalDurration = audio.nrFrames/48000;
mic_data = double(a);
fs = 48000;


		[b,a]=ellip(5,.2,80,[500]/(fs/2),'high');
		plot_data=mic_data./abs(max(mic_data));

        try
		[s,f,t]=fb_pretty_sonogram(filtfilt(b,a,mic_data./abs(max(mic_data))),fs,'low',2.5,'zeropad',0);
      
        
		minpt=1;
		maxpt=min(find(f>=10e3));

		%imwrite(flipdim(uint8(s(minpt:maxpt,:)),1),hot,fullfile(gif_dir,[file '.gif']),'gif');
		 catch
            disp('no audio... skipping spectrogram');
        end
        %save(fullfile(mat_dir,[file '.mat']),'audio','video','-v7.3');

        % clear the buffer
        
% Change audio data to proper frequency
audio = changeFrameRate(audio.data, 48000, 44150);
audio_ts = getAudTS(audio, 44150);
% Align video and audio data???
alignmentPoints = alignTimes(audio_ts, v_ts);


[Y, M, D, H, MN, S] = getDateVec(filenames{i}); 

% create serial time stamp
refTime = datenum([1903 12 31 20 00 00]);
timeStamp = datenum([str2num(Y), str2num(M), str2num(D), str2num(H), str2num(MN), str2num(S)]);
serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);
 
name = formFileName(birdName, serialNum, Y, M, D, H, MN, S);

aud_file_path = strcat(specTarget, '/', name, '.wav'); 
mov_file_path = strcat(specTarget, '/', name, '.mov');
csv_file_path = strcat(specTarget, '/', name, '.csv');

copyfile(filenames{i}, mov_file_path); % Put .mov file into new directory
copyfile(csvfilenames{i}, csv_file_path); % *should* put .csv file in new directory...

audiowrite(char(aud_file_path), audio, 44150);
        
clear video;
clear audio;
clear a_ts;
clear a;
clear v_ts;
clear v;
clear audio; 
end
% save all those indices!
fprintf(1,'\n');
