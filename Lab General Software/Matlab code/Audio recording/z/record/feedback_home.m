caller='Feedback';
birdname='test';
close_record_window;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANNELS
micro=1;         % -1 = off
speaker_input=-1;  % -1 = off
channels=[];  % [] = empty, electrode channels to be recorded
channels_o=[1];
scanrate=40000; 
scanrate_o=40000;
%scanrate=8000; scanrate_o=8000;
trigger_int=0.01; %duration of trigger signal in seconds
pre_trigger_int=1;  %pre-trigger in seconds
mid_trigger_int=0.8; % trigger during recording
post_trigger_int=1; %post-trigger in seconds
amp_scale=3000; % sets the sensitivity of the microphone
amp_drop=1;  % the factor below which the signal on the microphone must drop to stop triggering


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULTS
points_per_pixel=3; % the number of samples per pixel for the plotting routine
do_save_default=1;  
do_plot_default=1;
directed_default=1;   % 1=directed song,  0=undirected song
savepath='/Matlab-6.5/Data/';
pfeedback_default=0.0;
input_device='winsound'; input_hwaddress='0';
output_device='winsound'; output_hwaddress='0';

% NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE  NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise_int=.1; % minimum noise duration, in seconds
cutout_noise=[4000 6500]; % lower and upper limit of noise cutout, in Hz
cutout_mic=[4500 6000]; % lower and upper limit of micro cutout, in Hz

% GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY GOLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
golay_order=8;   golay_samples=2^golay_order;% duration in samples of golay pulse
intergolay_int=.025; % interval between golay pulses, in s
pre_golay_int=noise_int+.2; % trigger-time of golay pulse
wait_post_golay=0.05; % wait-time after golay before peeking microphone, as small as possible

% FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filt_noise_b,filt_noise_a]=ellip(6,2,60,cutout_noise/scanrate,'stop');
[filt_mic_b,filt_mic_a]=ellip(6,2,60,cutout_mic/scanrate);
 
define_dependent_variables;
new_noise_golay;
control_record_figure;
if do_plot
    control_record_figure_plot;
end
startup_engines;
