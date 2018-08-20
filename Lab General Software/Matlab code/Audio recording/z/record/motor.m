micro=1;         % -1 = off
loudspeaker=-1;  % -1 = off
channels=[];  % [] = empty, electrode channels to be recorded
do_save=1;   birdname='test';
%scanrate=44100;
scanrate=8000;
trigger_int=.01; % in seconds
pre_trigger_int=.5;  %pre-trigger in seconds
mid_trigger_int=0.5; % trigger during recording
post_trigger_int=1; %post-trigger in seconds
directed=1;   % 1=directed song,  0=undirected song
state_ai=0;  % 1=running    0=stop
amp_scale=120; % sets the sensitivity of the microphone
amp_drop=3;  % the factor below which the signal on the microphone must drop to stop triggering


caller='Motor';
dirname= datestr(date,29);  % directory for saving of data
thr_amp=3;% trigger_threshold (scale), dot not change this value
channels=sort(channels);
chans=[]; chans_col=[];  % define channel list and colors
if micro>-1
    chans=micro; chans_col='k'; chans_legend{1}='Microphone';
end 
if length(channels)>0
for i=1:length(chans)
    chans=[chans channels]; chans_col=[chans_col 'b']; 
    chans_legend{i+(micro>-1)}=['Channel ' num2str(channels(i))];
end
end
if loudspeaker>-1
    chans=[chans loudspeaker]; chans_col=[chans_col 'r'];
    chans_legend{length(chans_legend)+1}='Loudspeaker';
end
undirected=~directed;
max_amp=2*thr_amp; % maximum amplitude for display
axis_amp=2*thr_amp; % axis limits for display
curr_amp=.0;  % current_amplitude

pre_trigger_samples=ceil(pre_trigger_int*scanrate);
mid_trigger_samples=ceil(mid_trigger_int*scanrate);
post_trigger_samples=ceil(post_trigger_int*scanrate);
trigger_samples=ceil(trigger_int*scanrate);

close all;
ai = analoginput('winsound',0);
addchannel(ai, chans);
% Configure the analog input object.
set(ai, 'SampleRate', scanrate);


if do_save
    if ~(exist(birdname)==7)
        mkdir('.',birdname);
    end
    if ~(exist([birdname '/' dirname])==7)
        mkdir(['./' birdname],dirname);
    end
    if directed
        nameprefix='D';
    else
        nameprefix='U';
    end
    % look for existing files
   fnames=get_fnames(birdname,dirname,nameprefix,5000);
   if length(fnames)>0
        namehelp=char(fnames{end}(1:end-4));
        namelength=length([birdname nameprefix]);
        namehelp=NextName(namehelp);
  else
      namehelp=[birdname nameprefix '0001'];
    end
        filename=[birdname '/' dirname '/' namehelp];
    
    set(ai,'LogFileName',filename);
    set(ai,'LogToDiskMode','index');
    set(ai,'LoggingMode','Disk&Memory');
end 
% Configure the analog input object to trigger manually twice.

% Initialize callback parameters.  The TimerAction is initialized 
% after figure has been created.
set(ai,'StartFcn','t_start=clock;');
set(ai, 'SamplesPerTrigger',inf);
set(ai, 'TriggerRepeat', 0);
set(ai, 'TriggerType', 'manual');
set(ai, 'TriggerDelay', -pre_trigger_int);
set(ai, 'TriggerDelayUnits', 'seconds');
set(ai, 'TimerPeriod', 0.01);  

init_figure;
start(ai);
pause(trigger_int);
set(ai, 'TimerFcn', 'pre_trig');

%get(ai,'SamplesAvailable')
%get(ai,'SamplesAcquired')

%x = peekdata(ai, trigger_samples);

% Obtain the available time and data.
%[d,time] = getdata(ai, ai.SamplesPerTrigger);
