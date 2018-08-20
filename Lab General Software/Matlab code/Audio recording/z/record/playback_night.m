caller='Playback';
birdname='juvr6';
songtypeP=1;  % 1=BOS
is_forward=1;   % 1=forward,  0=reverse
close all;
bytes=8000;
path(path,'..');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANNELS
micro=8;         % -1 = off
speaker_input=-1;  % -1 = off
channels=[];  % [] = empty, electrode channels to be recorded
channels_o=[0];
scanrate=40000; 
scanrate_o=40000;

pre_bos_int=0; % time before BOS playbacks in seconds
post_bos_int=0; % time after BOS playbacks in seconds
playback_int=60;  % time between playbacks, in seconds
time_on=.85; % 0.85= 8:24pm
%time_on=.35; % 0.85= 8:24pm
time_off=.2; % 0.3=4:48 AM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULTS
points_per_pixel=3; % the number of samples per pixel for the plotting routine
loadpath='/Matlab-6.5/Data/';  % where BOS files are
savepath='/Matlab-6.5/Data/';  % where to save the data
input_device='nidaq'; input_hwaddress='1';
output_device='nidaq'; output_hwaddress='1';


load([loadpath birdname '/BOS_' birdname '-2003-02-27+noise.mat']);
sdata=lp;
if is_forward
bos_data=[zeros(1,pre_bos_int*scanrate) sdata zeros(1,post_bos_int*scanrate)]';
else
bos_data=[zeros(1,pre_bos_int*scanrate) sdata(end:-1:1) zeros(1,post_bos_int*scanrate)]';
end
bos_samples=length(bos_data);

daqreset
ao=analogoutput(output_device,output_hwaddress);
addchannel(ao,channels_o,'Loudspeaker');
set(ao,'SampleRate',scanrate_o);
set(ao,'TriggerType','Manual');
set(ao, 'RepeatOutput', 0);

putdata(ao,bos_data);
%fprintf('Press enter when ready for %s-%s playback to %s\n',songtypeP,caller,birdname);
%    fprintf('Next file: %s \n',namehelp);
    fprintf('Pause - hit a key\n');
pause
%%%%%%%%%%%%%%%%%%%%%%% Playback
ci=0;
while 1
if rem(now,1)>time_on | rem(now,1)<time_off
      start(ao); pause(1);
    while  strcmp(ao.Running,'Off')
    end
    trigger(ao);  
    ci=ci+1;fprintf('Currently: %d \n',ci);
      pause(playback_int-3-length(bos_data)/scanrate);
stop(ao); pause(1);
    putdata(ao,bos_data);
    pause(1);    
end
end
delete(ao);