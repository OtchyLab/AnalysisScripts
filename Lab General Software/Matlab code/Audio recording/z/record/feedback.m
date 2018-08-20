caller='Feedback';
birdname='juvenile_1'; %put the name of the file here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANNELS
micro=8;         % -1 = off
speaker_input=-1;  % -1 = off
channels=[];  % [] = empty, electrode channels to be recorded
channels_o=[0];
scanrate=30000; 
scanrate_o=30000;
%scanrate=8000; scanrate_o=8000;
peek_period=0.04; %duration of peekdata in seconds
rep_fcn_period=0.02; %period of rep_fcn execution in seconds
pre_trigger_int=1;  %pre-trigger in seconds
mid_trigger_int=0.8; % trigger during recording
post_trigger_int=1.2; %post-trigger in seconds
on_amp_scale=200; % sets the ON-sensitivity of the microphone
%on_amp_scale=200; % sets the ON-sensitivity of the microphone
off_amp_scale=150; % sets the OFF-sensitivity of the microphone

amp_drop=1;  % the factor below which the signal on the microphone must drop to stop triggering


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULTS
points_per_pixel=10; % the number of samples per pixel for the plotting routine
do_save_default=1;  
do_plot_default=1;
directed_default=1;   % 1=directed song,  0=undirected song
savepath='/Matlab-6.5/Data/';
pfeedback_default=0.0;
input_device='nidaq'; input_hwaddress='1';
output_device='nidaq'; output_hwaddress='1';

% Trigger windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
on_mic=[5200 5600]; % ON-window of microphone, in Hz
off_mic=[2000 4000]; % OFF-window of microphone, in Hz
do_off_trig=0;

% NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise_int=.1; % minimum noise duration, in seconds
cutout_noise=[4700 6100]; % lower and upper limit of noise cutout, in Hz

% GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
golay_order=8;   golay_samples=2^golay_order;% duration in samples of golay pulse
intergolay_int=.025; % interval between golay pulses, in s
pre_golay_int=noise_int+.4; % trigger-time of golay pulse
wait_post_golay=0.1; % wait-time after golay before peeking microphone, as small as possible


% FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filt_noise_b,filt_noise_a]=ellip(6,2,60,cutout_noise/scanrate*2,'stop');
[on_mic_b,on_mic_a]=ellip(6,2,60,on_mic/scanrate*2);
[off_mic_b,off_mic_a]=ellip(6,2,60,off_mic/scanrate*2);
 

define_dependent_variables;
new_noise_golay;
startup_engines;
control_record_figure;
if do_plot
    control_record_figure_plot;
end
